#!/usr/bin/env python3
"""Stage 12 — novelty verdict. Genome ANI (skani) is the strong signal when a reference
genome set is staged; otherwise the ITS-to-UNITE distance from Stage 08 is the fast
fallback. Thresholds: genome ANI < 95%  OR  best ITS < 98.5%  => candidate novel species,
reported with the GCPSR/polyphasic caveat (confirm by multi-locus + a genome ANI DB)."""
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
                ref = p[ci.get("Ref_name", 1)] if "Ref_name" in ci else (p[1] if len(p) > 1 else "?")
                if best is None or ani > best[0]:
                    best = (ani, ref)
    except Exception:
        return None
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--skani")
    ap.add_argument("--identify-json")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    skani = top_skani(a.skani)
    its_pid, species = None, None
    if a.identify_json and os.path.exists(a.identify_json):
        try:
            d = json.load(open(a.identify_json))
            its_pid = d.get("its_identity")
            species = d.get("species")
        except Exception:
            pass

    out = {"sample": a.sample, "stage": "novelty", "species": species,
           "genome_ani": None, "its_identity": its_pid}
    if skani:
        ani, ref = skani
        out["genome_ani"] = round(ani, 2); out["closest_reference"] = ref
        out["verdict"] = "candidate_novel_species" if ani < 95 else "known_species"
        out["basis"] = "genome_ani(skani)"
        out["note"] = ("ANI < 95% to closest reference genome — candidate novel species; confirm "
                       "by multi-locus GCPSR + polyphasic description") if ani < 95 else \
                      "matches a described species at species-level genome ANI"
    elif its_pid is not None:
        novel = its_pid < 98.5
        out["verdict"] = "candidate_novel_species" if novel else "known_species"
        out["basis"] = "ITS_distance(UNITE)"
        out["note"] = (f"best ITS identity {its_pid}% (< 98.5% species threshold) — candidate novel; "
                       "confirm with genome ANI + multi-locus GCPSR") if novel else \
                      f"ITS identity {its_pid}% is within the species range"
    else:
        out["verdict"] = "undetermined"
        out["basis"] = "none"
        out["note"] = "no ITS hit and no genome-ANI reference set staged"

    json.dump(out, open(a.out, "w"), indent=2)
    print(f"[novelty_call] {a.sample}: {out['verdict']} ({out['basis']})")


if __name__ == "__main__":
    main()
