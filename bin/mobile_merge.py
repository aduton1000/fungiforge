#!/usr/bin/env python3
"""Stage 10 merge (W2.8) — TE landscape (RepeatMasker table + RepeatModeler library), geNomad
(viruses / plasmids / proviruses on nuclear + mitochondrial contigs), the RVDB mycovirus / EVE
screen and the mitochondrial homing-endonuclease ORFs into one mobile.json."""
from __future__ import annotations
import argparse, csv, glob, json, os


def load(path):
    if path and os.path.exists(path):
        try:
            return json.load(open(path))
        except Exception:  # noqa: BLE001
            return {}
    return {}


def genomad_summary(gdir):
    """Counts from geNomad's *_summary/*_virus_summary.tsv and *_plasmid_summary.tsv (+ provirus flag)."""
    out = {"n_virus": 0, "n_plasmid": 0, "n_provirus": 0, "virus_taxa": {}, "ran": False}
    if not gdir or not os.path.isdir(gdir):
        return out
    for f in glob.glob(os.path.join(gdir, "**", "*_virus_summary.tsv"), recursive=True):
        out["ran"] = True
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out["n_virus"] += 1
                if (r.get("topology") or "").lower().startswith("provirus"):
                    out["n_provirus"] += 1
                tax = (r.get("taxonomy") or "").split(";")[-1] or "unclassified"
                out["virus_taxa"][tax] = out["virus_taxa"].get(tax, 0) + 1
    for f in glob.glob(os.path.join(gdir, "**", "*_plasmid_summary.tsv"), recursive=True):
        out["ran"] = True
        with open(f) as fh:
            out["n_plasmid"] += max(0, sum(1 for _ in fh) - 1)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True); ap.add_argument("--te"); ap.add_argument("--tbl"); ap.add_argument("--genomad")
    ap.add_argument("--mycovirus"); ap.add_argument("--heg"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    te, tbl, myco, heg = load(a.te), load(a.tbl), load(a.mycovirus), load(a.heg)
    gn = genomad_summary(a.genomad)
    out = {"sample": a.sample, "stage": "mobile",
           "te_percent": tbl.get("te_percent"), "te_landscape": tbl.get("summary"), "te_classes": tbl.get("classes"),
           "n_te_families": te.get("n_families", "NA"), "te_by_class": te.get("by_class", {}),
           "genomad": gn, "n_genomad_virus": gn["n_virus"] if gn["ran"] else None, "n_genomad_plasmid": gn["n_plasmid"] if gn["ran"] else None,
           "mycovirus": {k: v for k, v in myco.items() if k not in ("sample",)} if myco else None,
           "n_mycovirus": myco.get("n_mycovirus_eve_candidates") if myco else None,
           "mito_heg": {k: v for k, v in heg.items() if k not in ("sample",)} if heg else None,
           "n_mito_heg": heg.get("n_heg_orfs") if heg else None,
           "note": "mycovirus detection is DNA-only: RNA mycoviruses appear only as endogenous viral elements; retroelement-like RVDB hits are LTR retrotransposons"}
    json.dump(out, open(a.out, "w"), indent=2)
    print(f"[mobile_merge] {a.sample}: te_percent={out['te_percent']} families={out['n_te_families']} genomad_virus={out['n_genomad_virus']} "
          f"mycovirus_EVE={out['n_mycovirus']} mito_HEG={out['n_mito_heg']}")


if __name__ == "__main__":
    main()
