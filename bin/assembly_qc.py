#!/usr/bin/env python3
"""Stage 05 helper — assembly contiguity stats + BUSCO/compleasm completeness -> JSON.
Pure stdlib so it runs in any stage container. Sets a MIMAG-style qc_pass flag."""
from __future__ import annotations
import argparse, glob, json, os, re


def contig_lengths(fasta):
    lens, cur = [], 0
    if not fasta or not os.path.exists(fasta):
        return lens
    for line in open(fasta):
        if line.startswith(">"):
            if cur: lens.append(cur); cur = 0
        else:
            cur += len(line.strip())
    if cur: lens.append(cur)
    return lens


def n50(lens):
    if not lens: return 0
    s = sorted(lens, reverse=True); half = sum(s) / 2; acc = 0
    for L in s:
        acc += L
        if acc >= half: return L
    return s[-1]


def busco_complete(compleasm_dir, busco_dir):
    """Return % complete from compleasm summary or BUSCO short_summary, or None."""
    for p in glob.glob(os.path.join(compleasm_dir or "", "summary.txt")):
        t = open(p).read()
        m = re.search(r"S:\s*([0-9.]+)%", t)
        d = re.search(r"D:\s*([0-9.]+)%", t)
        if m:
            return round(float(m.group(1)) + (float(d.group(1)) if d else 0.0), 2)
    for p in glob.glob(os.path.join(busco_dir or "", "short_summary*.txt")):
        t = open(p).read()
        m = re.search(r"C:\s*([0-9.]+)%", t)
        if m:
            return float(m.group(1))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--nuclear", required=True)
    ap.add_argument("--lineage", default="fungi_odb10")
    ap.add_argument("--compleasm-dir", default="compleasm")
    ap.add_argument("--busco-dir", default="busco")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    lens = contig_lengths(a.nuclear)
    comp = busco_complete(a.compleasm_dir, a.busco_dir)
    total = sum(lens)
    out = {"sample": a.sample, "stage": "assembly_qc", "lineage": a.lineage,
           "assembly_len": total, "n_contigs": len(lens),
           "n50": n50(lens), "largest_contig": max(lens) if lens else 0,
           "busco_complete": comp,
           "qc_pass": bool(total > 0 and (comp is None or comp >= 80.0))}
    json.dump(out, open(a.out, "w"), indent=2)
    print(f"[assembly_qc] {a.sample}: {total/1e6:.2f} Mb, {len(lens)} contigs, "
          f"N50 {out['n50']}, BUSCO {comp if comp is not None else 'NA'}%")


if __name__ == "__main__":
    main()
