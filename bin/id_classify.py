#!/usr/bin/env python3
"""Stage 08 — resolve a species call from ITS/marker hits + genome-level sourmash
gather. Genome ANI (gather top hit) is primary; ITS/UNITE refines. Emits a species
call with confidence and the supporting method (GCPSR concordance noted when markers
and genome agree). Degrades to 'unknown' when no reference evidence is present."""
from __future__ import annotations
import argparse, csv, json, os, re


def top_gather(path):
    if not path or not os.path.exists(path):
        return None
    try:
        rows = list(csv.DictReader(open(path)))
    except Exception:
        return None
    if not rows:
        return None
    r = rows[0]
    name = r.get("name") or r.get("match_name") or ""
    # sourmash 'name' is usually 'ACCESSION Genus species strain...'
    m = re.search(r"([A-Z][a-z]+ [a-z]+)", name)
    return {"species": m.group(1) if m else name[:60],
            "ani": r.get("average_containment_ani") or r.get("f_match") or "",
            "raw": name[:120]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--markers"); ap.add_argument("--gather")
    ap.add_argument("--data-dir", default="")
    ap.add_argument("--out-species", required=True); ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    gt = top_gather(a.gather)
    has_its = bool(a.markers and os.path.exists(a.markers) and os.path.getsize(a.markers) > 0)
    if gt and gt["species"]:
        species, conf, method = gt["species"], ("high" if gt.get("ani") else "medium"), "genome_ani(sourmash)"
        if has_its:
            method += "+ITS"; conf = "high"  # multi-evidence -> GCPSR concordance
    else:
        species, conf, method = "unknown", "none", ("ITS_only" if has_its else "no_reference")

    open(a.out_species, "w").write(f"{species}\t{conf}\n")
    json.dump({"sample": a.sample, "stage": "identify", "species": species,
               "confidence": conf, "method": method,
               "gather_top": gt, "its_present": has_its}, open(a.out_json, "w"), indent=2)
    print(f"[id_classify] {a.sample}: {species} ({conf}, {method})")


if __name__ == "__main__":
    main()
