#!/usr/bin/env python3
"""Stage 10 merge — combine the TE summary with the geNomad mycovirus/EVE scan into
one mobile-elements result.json (DNA-only caveat noted for mycoviruses)."""
from __future__ import annotations
import argparse, glob, json, os


def genomad_counts(gdir):
    n_vir = 0
    for f in glob.glob(os.path.join(gdir or "", "**", "*virus_summary.tsv"), recursive=True):
        try:
            n_vir += max(0, sum(1 for _ in open(f)) - 1)
        except Exception:
            pass
    return n_vir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--te"); ap.add_argument("--genomad")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    te = {}
    if a.te and os.path.exists(a.te):
        try: te = json.load(open(a.te))
        except Exception: te = {}
    out = {"sample": a.sample, "stage": "mobile",
           "n_te_families": te.get("n_families", "NA"), "te_by_class": te.get("by_class", {}),
           "n_mycovirus": genomad_counts(a.genomad),
           "note": "mycovirus detection is DNA-only (RNA mycoviruses not captured by WGS)"}
    json.dump(out, open(a.out, "w"), indent=2)
    print(f"[mobile_merge] {a.sample}: {out['n_te_families']} TE families, {out['n_mycovirus']} viral contigs")


if __name__ == "__main__":
    main()
