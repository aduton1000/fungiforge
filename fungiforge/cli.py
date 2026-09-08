"""FungiForge command-line front door (stdlib-only, zero heavy deps).

Subcommands:
  run          wrap `nextflow run main.nf` with sane defaults
  samplesheet  build a samplesheet CSV from directories of ONT and/or Illumina reads
               (ONT-only, Illumina-only, or hybrid — per sample, whichever are found)
  fetch-refs   print the reference-database fetch command (bin/fetch_references.sh)
  version
"""
from __future__ import annotations
import argparse
import glob
import os
import subprocess
import sys
from . import __version__

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def cmd_run(a):
    cmd = ["nextflow", "run", os.path.join(REPO, "main.nf"), "-profile", a.profile]
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
        s = sid(p)
        for tag in ("_R1", "_R2", "_1", "_2"):
            s = s.replace(tag, "")
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
    print(f'bash {os.path.join(REPO, "bin", "fetch_references.sh")}')


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

    v = sub.add_parser("version"); v.set_defaults(func=cmd_version)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
