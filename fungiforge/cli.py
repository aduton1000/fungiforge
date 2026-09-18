"""FungiForge command-line front door (stdlib-only, zero heavy deps).

Subcommands:
  run          wrap `nextflow run main.nf` with sane defaults
  samplesheet  build a samplesheet CSV from directories of ONT and/or Illumina reads
               (ONT-only, Illumina-only, or hybrid — per sample, whichever are found)
  fetch-refs   print the reference-database fetch command (bin/fetch_references.sh)
  validate     compare a finished run with an expected-results file (bin/validate_run.py)
  check        validate a samplesheet before a run (bin/validate_samplesheet.py)
  check-master validate a master table against the schema (bin/validate_master.py)
  version
"""
from __future__ import annotations
import argparse
import re
import glob
import os
import shutil
import subprocess
import sys
from . import __version__


def _repo_root():
    """Locate the pipeline checkout (main.nf, bin/, test/expected/).

    Deriving it from __file__ only works when the CLI runs from a source tree. A site
    install pip-installs the package into cli-env/lib/.../site-packages, where the
    parent directory holds no bin/, and every subcommand that shells out to a helper
    died on a path that cannot exist. FUNGIFORGE_HOME is set by the site environment
    (share/fungiforge-env.sh.example) and names the checkout, so prefer it.
    """
    for cand in (os.environ.get("FUNGIFORGE_HOME"), os.environ.get("FUNGIFORGE_ROOT", "") and
                 os.path.join(os.environ["FUNGIFORGE_ROOT"], "repo"),
                 os.path.dirname(os.path.dirname(os.path.abspath(__file__)))):
        if cand and os.path.isfile(os.path.join(cand, "main.nf")):
            return os.path.abspath(cand)
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


REPO = _repo_root()


def helper(name):
    """Absolute path to a bin/ helper script, or exit with an actionable message.

    Falls back to PATH: the site environment puts $FUNGIFORGE_HOME/bin on it, so an
    installed CLI outside a checkout still finds the scripts.
    """
    local = os.path.join(REPO, "bin", name)
    if os.path.isfile(local):
        return local
    found = shutil.which(name)
    if found:
        return found
    sys.exit("[fungiforge] cannot find bin/%s. Looked in %s and on PATH. Point FUNGIFORGE_HOME "
             "at the pipeline checkout (the directory holding main.nf), or run from it." % (name, REPO))


def cmd_run(a):
    entry = os.path.join(REPO, "main.nf")
    if not os.path.isfile(entry):
        sys.exit("[fungiforge] cannot find main.nf under %s. Point FUNGIFORGE_HOME at the "
                 "pipeline checkout." % REPO)
    cmd = ["nextflow", "run", entry, "-profile", a.profile]
    if a.samplesheet: cmd += ["--samplesheet", a.samplesheet]
    if a.data_dir:    cmd += ["--data_dir", a.data_dir]
    if a.outdir:      cmd += ["--outdir", a.outdir]
    if a.resume:      cmd += ["-resume"]
    cmd += a.extra
    print("[fungiforge]", " ".join(cmd))
    return subprocess.call(cmd)


def cmd_samplesheet(a):
    """Build a samplesheet from ONT and/or Illumina reads, keyed by sample id.

    A sample may have ONT only, Illumina only, or both (hybrid) — every sample id
    seen in any of the three inputs gets a row; missing columns are left blank and
    the pipeline routes each isolate by its inputs (see main.nf meta.assembly_mode).
    """
    def sid(p):
        b = os.path.basename(p)
        for suf in (".ont.fastq.gz", ".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if b.endswith(suf): return b[: -len(suf)]
        return os.path.splitext(b)[0]

    def ilmn_sid(p):
        # Strip the read-pair designation without touching the rest of the name:
        #   bcl2fastq/BCL Convert  NAME_S12_R1_001  -> NAME   (also NAME_S12_L001_R1_001)
        #   SRA / generic          NAME_R1, NAME_1, NAME.R1 -> NAME
        s = sid(p)
        s = re.sub(r"_S\d+(_L\d{3})?_R[12]_\d{3}$", "", s)
        s = re.sub(r"[._]R?[12]$", "", s)
        return s

    ont = {sid(p): p for pat in (a.ont or []) for p in sorted(glob.glob(pat))}
    r1 = {ilmn_sid(p): p for pat in (a.illumina_r1 or []) for p in sorted(glob.glob(pat))}
    r2 = {ilmn_sid(p): p for pat in (a.illumina_r2 or []) for p in sorted(glob.glob(pat))}
    if not (ont or r1 or r2):
        sys.exit("[fungiforge] no reads matched --ont / --illumina-r1 / --illumina-r2 globs")

    samples = sorted(set(ont) | set(r1) | set(r2))
    rows = ["sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season"]
    n_ont = n_ilmn = n_hyb = 0
    dropped = []
    for s in samples:
        o, a1, a2 = ont.get(s, ""), r1.get(s, ""), r2.get(s, "")
        has_ont, has_pair = bool(o), bool(a1 and a2)
        if not has_ont and not has_pair:
            dropped.append(s)  # e.g. an unpaired R1 with no ONT — cannot assemble
            continue
        if has_ont and has_pair: n_hyb += 1
        elif has_ont:            n_ont += 1
        else:                    n_ilmn += 1
        rows.append(f"{s},{o},{a1},{a2},{a.compartment},{a.facility},{a.season}")
    out = a.out or "samplesheet.csv"
    open(out, "w").write("\n".join(rows) + "\n")
    print(f"[fungiforge] {len(rows)-1} samples -> {out} "
          f"(ONT-only={n_ont}, Illumina-only={n_ilmn}, hybrid={n_hyb})")
    if dropped:
        print(f"[fungiforge] WARNING: skipped {len(dropped)} sample(s) with neither ONT nor "
              f"paired Illumina: {', '.join(dropped)}")


def cmd_fetch_refs(a):
    print(f'export FUNGIFORGE_DB="{a.data_dir}"')
    print(f'bash {helper("fetch_references.sh")}')


def cmd_validate(a):
    """Compare results/ with test/expected/<sample>.json (exit 1 on any failed check)."""
    exp = a.expected or os.path.join(REPO, "test", "expected", f"{a.sample}.json")
    cmd = [sys.executable, helper("validate_run.py"), "--results", a.results, "--expected", exp]
    if a.sample: cmd += ["--sample", a.sample]
    if a.out:    cmd += ["--out", a.out]
    return subprocess.call(cmd)


def cmd_check(a):
    """Validate a samplesheet before launching (exit 1 on any error)."""
    cmd = [sys.executable, helper("validate_samplesheet.py"), "--samplesheet", a.samplesheet]
    if a.json: cmd += ["--json", a.json]
    if a.strict: cmd += ["--strict"]
    if a.no_check_files: cmd += ["--no-check-files"]
    return subprocess.call(cmd)


def cmd_check_master(a):
    """Validate a master table against fungiforge/resources/master_schema.json."""
    cmd = [sys.executable, helper("validate_master.py"), "--master", a.master]
    if a.json: cmd += ["--json", a.json]
    if a.strict: cmd += ["--strict"]
    return subprocess.call(cmd)


def cmd_version(a):
    print(f"FungiForge v{__version__}")


def build_parser():
    p = argparse.ArgumentParser(prog="fungiforge", description="FungiForge — fungal ONT / Illumina / hybrid genomics")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run the Nextflow pipeline")
    r.add_argument("--samplesheet"); r.add_argument("--data_dir"); r.add_argument("--outdir")
    r.add_argument("-profile", "--profile", dest="profile", default="local,docker")
    r.add_argument("-resume", "--resume", action="store_true")
    r.add_argument("extra", nargs=argparse.REMAINDER, help="extra args passed through to nextflow")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("samplesheet", help="build a samplesheet from read directories (ONT and/or Illumina)")
    s.add_argument("--ont", nargs="*", help="ONT fastq files/globs (optional if Illumina given)")
    s.add_argument("--illumina-r1", nargs="*", help="Illumina R1 fastq files/globs")
    s.add_argument("--illumina-r2", nargs="*", help="Illumina R2 fastq files/globs")
    s.add_argument("--compartment", default="NA"); s.add_argument("--facility", default="NA"); s.add_argument("--season", default="NA")
    s.add_argument("-o", "--out")
    s.set_defaults(func=cmd_samplesheet)

    f = sub.add_parser("fetch-refs", help="print the reference-DB fetch command")
    f.add_argument("--data_dir", required=True)
    f.set_defaults(func=cmd_fetch_refs)

    va = sub.add_parser("validate", help="compare a finished run with an expected-results file")
    va.add_argument("--results", required=True, help="the run's --outdir")
    va.add_argument("--sample", default=None, help="sample id (default: from the expected file)")
    va.add_argument("--expected", default=None, help="expected JSON (default: test/expected/<sample>.json)")
    va.add_argument("-o", "--out", default=None, help="write the check table as TSV")
    va.set_defaults(func=cmd_validate)

    c = sub.add_parser("check", help="validate a samplesheet before a run")
    c.add_argument("--samplesheet", required=True)
    c.add_argument("--json", default=None, help="write the report as JSON")
    c.add_argument("--strict", action="store_true", help="treat warnings as errors")
    c.add_argument("--no-check-files", action="store_true", help="do not open the read files (headers and ids only)")
    c.set_defaults(func=cmd_check)

    cm = sub.add_parser("check-master", help="validate a master table against the schema")
    cm.add_argument("--master", required=True, help="master_fungi.tsv or a per-isolate *.master.tsv")
    cm.add_argument("--json", default=None); cm.add_argument("--strict", action="store_true")
    cm.set_defaults(func=cmd_check_master)

    v = sub.add_parser("version"); v.set_defaults(func=cmd_version)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
