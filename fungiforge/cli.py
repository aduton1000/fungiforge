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
import textwrap
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


class _Help(argparse.RawDescriptionHelpFormatter):
    """Keep hand-written sections verbatim, and drop argparse's `{a,b,c}` metavar line.

    That brace list is the single ugliest thing in a default argparse CLI: it repeats every
    command name in a line nobody reads, wraps badly once there are more than four, and pushes
    the real descriptions down. The subcommand descriptions below it say the same thing properly.
    """

    def __init__(self, prog):
        super().__init__(prog, max_help_position=30, width=88)

    def add_argument(self, action):
        super().add_argument(action)
        # argparse sizes the description column from the longest invocation it has seen, which
        # leaves the longest command name one character short of fitting and wraps it onto its
        # own line. Reserve the room up front so the column is straight.
        self._action_max_length = max(self._action_max_length, 18)

    def _format_action(self, action):
        text = super()._format_action(action)
        if action.nargs == argparse.PARSER:
            text = "".join(text.split("\n", 1)[1:])
        return text


def _d(text):
    """Wrap a description to the help width.

    The formatter keeps hand-written text verbatim so the examples below stay laid out, which
    means descriptions have to arrive already wrapped or they run off the terminal.
    """
    return textwrap.fill(" ".join(text.split()), 86)


EPILOG = """\
examples
  fungiforge samplesheet --ont 'reads/*.fastq.gz' -o samples.csv
  fungiforge check --samplesheet samples.csv
  fungiforge run --samplesheet samples.csv --outdir results -resume
  fungiforge validate --results results --sample AfumCEA10

on a cluster
  Use the site launcher rather than `fungiforge run`: it resolves the profile,
  databases, container cache and licensed paths from the site environment.
      fungiforge-run --samplesheet samples.csv --outdir results -resume

documentation
  docs/manual/fungiforge_manual.md       full manual, stage by stage
  docs/hpc_deployment.md                 install and database staging
  fungiforge <command> --help            options for one command
"""


def build_parser():
    p = argparse.ArgumentParser(
        prog="fungiforge", formatter_class=_Help, epilog=EPILOG,
        usage="fungiforge <command> [options]",
        description="FungiForge %s — whole-genome analysis of fungal isolates from ONT,\n"
                    "Illumina or hybrid reads: assembly, identification, antifungal\n"
                    "resistance, gene clusters and cohort phylogenomics." % __version__)
    p.add_argument("-V", "--version", action="version", version="FungiForge %s" % __version__,
                   help="print the version and exit")
    sub = p.add_subparsers(dest="cmd", required=True, title="commands", metavar="<command>",
                           prog="fungiforge")

    r = sub.add_parser("run", formatter_class=_Help, help="launch the Nextflow pipeline",
                       description=_d("Run the pipeline. Every unrecognised argument is passed "
                                   "through to `nextflow run` unchanged."),
                       epilog="example\n  fungiforge run --samplesheet samples.csv "
                              "--outdir results -resume\n")
    r.add_argument("--samplesheet", metavar="FILE", help="sample sheet CSV (see the manual)")
    r.add_argument("--data_dir", metavar="DIR", help="reference-database root")
    r.add_argument("--outdir", metavar="DIR", help="where results are written")
    r.add_argument("-profile", "--profile", dest="profile", default="local,docker", metavar="P",
                   help="Nextflow profile(s) (default: %(default)s)")
    r.add_argument("-resume", "--resume", action="store_true", help="reuse cached tasks")
    r.add_argument("extra", nargs=argparse.REMAINDER, metavar="...",
                   help="further arguments passed straight to nextflow")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("samplesheet", formatter_class=_Help,
                       help="build a sample sheet from directories of reads",
                       description=_d("Build a sample sheet from read files or globs. Isolate ids "
                                   "are taken from the file names; ONT-only, Illumina-only and "
                                   "hybrid rows may be mixed."),
                       epilog="example\n  fungiforge samplesheet --ont 'ont/*.fastq.gz' \\\n"
                              "      --illumina-r1 'ill/*_R1.fastq.gz' --illumina-r2 'ill/*_R2.fastq.gz' \\\n"
                              "      --compartment AIR --facility Abattoir_1 -o samples.csv\n")
    s.add_argument("--ont", nargs="*", metavar="GLOB", help="ONT FASTQ files or globs")
    s.add_argument("--illumina-r1", nargs="*", metavar="GLOB", help="Illumina R1 files or globs")
    s.add_argument("--illumina-r2", nargs="*", metavar="GLOB", help="Illumina R2 files or globs")
    s.add_argument("--compartment", default="NA", metavar="S", help="metadata column (default: %(default)s)")
    s.add_argument("--facility", default="NA", metavar="S", help="metadata column (default: %(default)s)")
    s.add_argument("--season", default="NA", metavar="S", help="metadata column (default: %(default)s)")
    s.add_argument("-o", "--out", metavar="FILE", help="write here instead of standard output")
    s.set_defaults(func=cmd_samplesheet)

    c = sub.add_parser("check", formatter_class=_Help,
                       help="validate a sample sheet before a run",
                       description=_d("Check a sample sheet before launching: duplicate or malformed "
                                   "ids, unpaired reads, missing or empty files, metadata typos. "
                                   "Catches at the sheet what would otherwise fail hours in."),
                       epilog="example\n  fungiforge check --samplesheet samples.csv --strict\n")
    c.add_argument("--samplesheet", required=True, metavar="FILE", help="the sheet to check")
    c.add_argument("--json", default=None, metavar="FILE", help="write the report as JSON")
    c.add_argument("--strict", action="store_true", help="treat warnings as errors")
    c.add_argument("--no-check-files", action="store_true",
                   help="check ids and pairing only, do not open the read files")
    c.set_defaults(func=cmd_check)

    va = sub.add_parser("validate", formatter_class=_Help,
                        help="compare a finished run against expected results",
                        description=_d("Compare a completed run with a recorded expectation "
                                    "(assembly size, contiguity, completeness, species, "
                                    "resistance, cluster count) and report each check."),
                        epilog="example\n  fungiforge validate --results results --sample AfumCEA10\n")
    va.add_argument("--results", required=True, metavar="DIR", help="the run's --outdir")
    va.add_argument("--sample", default=None, metavar="ID", help="isolate id (default: from the expected file)")
    va.add_argument("--expected", default=None, metavar="FILE",
                    help="expectation JSON (default: test/expected/<sample>.json)")
    va.add_argument("-o", "--out", default=None, metavar="FILE", help="write the check table as TSV")
    va.set_defaults(func=cmd_validate)

    cm = sub.add_parser("check-master", formatter_class=_Help,
                        help="validate a master table against the schema",
                        description=_d("Validate a master table against the published schema: column "
                                    "presence and order, types, enumerations, duplicate isolates "
                                    "and ragged rows."),
                        epilog="example\n  fungiforge check-master --master results/04_summary/master_fungi.tsv\n")
    cm.add_argument("--master", required=True, metavar="FILE",
                    help="master_fungi.tsv or a per-isolate *.master.tsv")
    cm.add_argument("--json", default=None, metavar="FILE", help="write the report as JSON")
    cm.add_argument("--strict", action="store_true", help="treat warnings as errors")
    cm.set_defaults(func=cmd_check_master)

    f = sub.add_parser("fetch-refs", formatter_class=_Help,
                       help="print the reference-database fetch command",
                       description=_d("Print the command that stages the reference databases. It is "
                                   "printed rather than run: it takes hours and belongs in a "
                                   "session you can detach from."),
                       epilog="example\n  fungiforge fetch-refs --data_dir /data/fungiforge\n")
    f.add_argument("--data_dir", required=True, metavar="DIR", help="where the databases go")
    f.set_defaults(func=cmd_fetch_refs)

    v = sub.add_parser("version", formatter_class=_Help, help="print the version",
                       description=_d("Print the FungiForge version."))
    v.set_defaults(func=cmd_version)
    return p


def main(argv=None):
    parser = build_parser()
    given = sys.argv[1:] if argv is None else list(argv)
    if not given:
        # A bare `fungiforge` should teach, not scold. argparse's "the following arguments are
        # required" tells a first-time user nothing about what the tool does.
        parser.print_help()
        return 2
    args = parser.parse_args(given)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
