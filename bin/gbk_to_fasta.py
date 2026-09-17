#!/usr/bin/env python3
"""Write the sequences of a GenBank file as FASTA, record ids as names (W2.4).

Funannotate renames and sorts contigs, so its GenBank coordinates do not map onto the
assembly FASTA the pipeline holds; reads for the read-level resistance check are therefore
mapped to the GenBank records themselves, and this writes that reference.
  gbk_to_fasta.py --gbk isolate.gbk --out records.fa
"""
from __future__ import annotations
import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gbk", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    from Bio import SeqIO
    n = 0
    with open(a.out, "w") as fh:
        for rec in SeqIO.parse(a.gbk, "genbank"):
            s = str(rec.seq)
            if not s or set(s.upper()) <= {"N"}:
                continue
            fh.write(f">{rec.id}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")
            n += 1
    print(f"[gbk_to_fasta] {n} records -> {a.out}")


if __name__ == "__main__":
    main()
