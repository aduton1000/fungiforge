"""FungiForge command-line front door (stdlib-only, zero heavy deps).

Subcommands:
  run          wrap `nextflow run main.nf` with sane defaults
  samplesheet  build a samplesheet CSV from directories of ONT (+ Illumina) reads
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
    """Pair ONT fastqs (required) with Illumina R1/R2 (optional) by sample id."""
    def sid(p):
        b = os.path.basename(p)
        for suf in (".ont.fastq.gz", ".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if b.endswith(suf): return b[: -len(suf)]
        return os.path.splitext(b)[0]

    ont = {sid(p): p for pat in a.ont for p in sorted(glob.glob(pat))}
    r1 = {sid(p).replace("_R1", "").replace("_1", ""): p for pat in (a.illumina_r1 or []) for p in sorted(glob.glob(pat))}
    r2 = {sid(p).replace("_R2", "").replace("_2", ""): p for pat in (a.illumina_r2 or []) for p in sorted(glob.glob(pat))}
    rows = ["sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season"]
    for s, p in ont.items():
        rows.append(f"{s},{p},{r1.get(s,'')},{r2.get(s,'')},{a.compartment},{a.facility},{a.season}")
    out = a.out or "samplesheet.csv"
    open(out, "w").write("\n".join(rows) + "\n")
    print(f"[fungiforge] {len(ont)} samples -> {out}")


def cmd_fetch_refs(a):
    print(f'export FUNGIFORGE_DB="{a.data_dir}"')
    print(f'bash {os.path.join(REPO, "bin", "fetch_references.sh")}')


def cmd_version(a):
    print(f"FungiForge v{__version__}")


def build_parser():
    p = argparse.ArgumentParser(prog="fungiforge", description="FungiForge — fungal ONT/hybrid genomics")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run the Nextflow pipeline")
    r.add_argument("--samplesheet"); r.add_argument("--data_dir"); r.add_argument("--outdir")
    r.add_argument("-profile", "--profile", dest="profile", default="local,docker")
    r.add_argument("-resume", "--resume", action="store_true")
    r.add_argument("extra", nargs=argparse.REMAINDER, help="extra args passed through to nextflow")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("samplesheet", help="build a samplesheet from read directories")
    s.add_argument("--ont", nargs="+", required=True, help="ONT fastq files/globs")
    s.add_argument("--illumina-r1", nargs="*"); s.add_argument("--illumina-r2", nargs="*")
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
