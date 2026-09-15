#!/usr/bin/env python3
"""Extract the rRNA-operon region(s) from a genome assembly using barrnap coordinates,
so ITSx sees SHORT sequences (a few kb) instead of whole chromosomes — ITSx/HMMER
aborts on sequences >100 kb ("over comparison pipeline limit"). We take a window around
every rRNA feature (18S / 5.8S / 28S) and merge overlaps; the ITS1-5.8S-ITS2 barcode
lives between the SSU and LSU, so a padded window captures it in full."""
from __future__ import annotations
import argparse


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name: seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--pad", type=int, default=3000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    contigs = read_fasta(a.assembly)
    # collect rRNA intervals per contig from the barrnap GFF3
    intervals = {}
    try:
        for line in open(a.gff):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[2] not in ("rRNA", "gene"):
                continue
            c, s, e = f[0], int(f[3]), int(f[4])
            intervals.setdefault(c, []).append((s, e))
    except Exception:
        pass

    n = 0
    with open(a.out, "w") as out:
        for c, ivs in intervals.items():
            seq = contigs.get(c, "")
            if not seq:
                continue
            # pad and merge overlapping windows
            wins = sorted((max(0, s - a.pad), min(len(seq), e + a.pad)) for s, e in ivs)
            merged = []
            for w in wins:
                if merged and w[0] <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], w[1]))
                else:
                    merged.append(list(w))
            for i, (s, e) in enumerate(merged):
                out.write(f">{c}_rRNAregion_{i}_{s}-{e}\n{seq[s:e]}\n")
                n += 1
    print(f"[extract_rrna_region] wrote {n} rRNA-operon window(s) -> {a.out}")


if __name__ == "__main__":
    main()
