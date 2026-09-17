#!/usr/bin/env python3
"""Collect an isolate's single-copy BUSCO proteins for cohort phylogenomics (W3.1, stage 05).

  busco_singlecopy.py --sample ID [--busco-dir busco] [--compleasm-dir compleasm] --out ID.busco_sc.faa

BUSCO layout:      <busco-dir>/run_*/busco_sequences/single_copy_busco_sequences/<id>.faa
compleasm layout:  <compleasm-dir>/<lineage>/full_table.tsv (Status == Single) + translated_protein.fasta
Output headers: ><busco_id>|<sample>. Empty file (and a note on stderr) when neither layout is found."""
from __future__ import annotations
import argparse, glob, os, re, sys


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def from_busco(busco_dir):
    out = {}
    for f in sorted(glob.glob(os.path.join(busco_dir or "", "run_*", "busco_sequences", "single_copy_busco_sequences", "*.faa"))):
        bid = os.path.basename(f)[:-4]
        s = read_fasta(f)
        if s:
            out[bid] = next(iter(s.values()))
    return out


def from_compleasm(comp_dir):
    out = {}
    for table in sorted(glob.glob(os.path.join(comp_dir or "", "*", "full_table.tsv"))):
        lineage_dir = os.path.dirname(table)
        prot_file = os.path.join(lineage_dir, "translated_protein.fasta")
        if not os.path.exists(prot_file):
            continue
        prots = read_fasta(prot_file)
        single = set()
        with open(table) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            ci = {c.strip("# "): i for i, c in enumerate(header)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > 1 and f[ci.get("Status", 1)] == "Single":
                    single.add(f[ci.get("Gene", 0)])
        for name, seq in prots.items():
            bid = re.sub(r"_\d+$", "", name)
            if bid in single and bid not in out:
                out[bid] = seq
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--busco-dir", default="busco"); ap.add_argument("--compleasm-dir", default="compleasm")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    seqs = from_compleasm(a.compleasm_dir) or from_busco(a.busco_dir)
    with open(a.out, "w") as fh:
        for bid in sorted(seqs):
            fh.write(f">{bid}|{a.sample}\n{seqs[bid].rstrip('*')}\n")
    if not seqs:
        sys.stderr.write(f"[busco_singlecopy] {a.sample}: no single-copy BUSCO sequences found (compleasm or BUSCO layout)\n")
    print(f"[busco_singlecopy] {a.sample}: {len(seqs)} single-copy BUSCO proteins -> {a.out}")


if __name__ == "__main__":
    main()
