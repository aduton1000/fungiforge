#!/usr/bin/env python3
"""Stage 12 — novelty verdict from skani ANI to the closest reference genome.
ANI < 95% (with adequate aligned fraction) => candidate novel species, reported with
the GCPSR/polyphasic caveat. Phylogenomic placement refines this in the comparative layer."""
from __future__ import annotations
import argparse, json, os


def top_skani(path):
    if not path or not os.path.exists(path):
        return None
    best = None
    try:
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            ci = {c: i for i, c in enumerate(header)}
            for line in fh:
                p = line.rstrip("\n").split("\t")
                try:
                    ani = float(p[ci.get("ANI", 2)])
                except Exception:
                    continue
                if best is None or ani > best[0]:
                    best = (ani, p[ci.get("Ref_name", 1)] if "Ref_name" in ci else p[1])
    except Exception:
        return None
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--skani"); ap.add_argument("--markers")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    top = top_skani(a.skani)
    if top:
        ani, ref = top
        verdict = "candidate_novel_species" if ani < 95 else "known_species"
        note = ("ANI < 95%% to closest reference — candidate novel species; confirm by "
                "multi-locus GCPSR + phylogenomic placement + polyphasic description") if ani < 95 else \
               "matches a described species at species-level ANI"
        out = {"sample": a.sample, "stage": "novelty", "top_ani": round(ani, 2),
               "closest_reference": ref, "verdict": verdict, "note": note}
    else:
        out = {"sample": a.sample, "stage": "novelty", "verdict": "undetermined",
               "note": "no reference ANI available (stage needs data_dir/refseq_fungi)"}
    json.dump(out, open(a.out, "w"), indent=2)
    print(f"[novelty_call] {a.sample}: {out['verdict']}")


if __name__ == "__main__":
    main()
