#!/usr/bin/env python3
"""Stage 08 — resolve a species call. Primary signal: ITS classified against UNITE
(vsearch blast6) — the standard fungal barcode, fast and emulation-friendly. Secondary:
genome-level sourmash gather (optional; slow under x86 emulation). GCPSR concordance is
reported when ITS and genome agree. Also emits the top ITS identity for novelty (Stage 12).

ITS species thresholds (guidance): >=98.5% -> species-level; 94-98.5% -> genus/uncertain;
<94% -> candidate novel (confirm with multi-locus + genome ANI)."""
from __future__ import annotations
import argparse, csv, json, os, re


def parse_unite_b6(path):
    """Best UNITE hit from a vsearch --blast6out file -> (species, genus, pident, sh)."""
    if not path or not os.path.exists(path):
        return None
    best = None
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            pid = float(f[2])
        except ValueError:
            continue
        target = f[1]
        if best is None or pid > best[0]:
            best = (pid, target)
    if not best:
        return None
    pid, target = best
    sp = re.search(r"s__([A-Za-z0-9_.\-]+)", target)
    ge = re.search(r"g__([A-Za-z0-9_.\-]+)", target)
    sh = re.search(r"(SH\d+\.\d+FU)", target)
    species = sp.group(1).replace("_", " ") if sp else None
    if species and re.search(r"_sp$|Incertae|unidentified|uncultured", species, re.I):
        species = None
    return {"pident": round(pid, 2),
            "species": species,
            "genus": ge.group(1) if ge else None,
            "sh": sh.group(1) if sh else None}


def top_gather(path):
    if not path or not os.path.exists(path):
        return None
    try:
        rows = list(csv.DictReader(open(path)))
    except Exception:
        return None
    if not rows:
        return None
    name = rows[0].get("name") or rows[0].get("match_name") or ""
    m = re.search(r"([A-Z][a-z]+ [a-z]+)", name)
    return {"species": m.group(1) if m else name[:60], "raw": name[:120]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--its")                  # ITS fasta (presence flag)
    ap.add_argument("--unite-b6")             # vsearch blast6 vs UNITE
    ap.add_argument("--gather")               # optional sourmash gather csv
    ap.add_argument("--out-species", required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    its_present = bool(a.its and os.path.exists(a.its) and os.path.getsize(a.its) > 0)
    unite = parse_unite_b6(a.unite_b6)
    gather = top_gather(a.gather)

    species, conf, method, its_pid = "unknown", "none", "no_reference", None
    if unite:
        its_pid = unite["pident"]
        if unite["species"] and its_pid >= 98.5:
            species, conf, method = unite["species"], "high", "ITS/UNITE"
        elif unite["species"] and its_pid >= 94:
            species, conf, method = unite["species"], "medium", "ITS/UNITE(below-species-threshold)"
        elif unite["genus"]:
            species, conf, method = unite["genus"] + " sp.", "low", "ITS/UNITE(genus-only)"
    # genome-level concordance / fallback
    if gather and gather["species"]:
        if species != "unknown" and gather["species"].split()[0] == species.split()[0]:
            method += "+genome_ani(concordant/GCPSR)"; conf = "high"
        elif species == "unknown":
            species, conf, method = gather["species"], "medium", "genome_ani(sourmash)"

    open(a.out_species, "w").write(f"{species}\t{conf}\n")
    json.dump({"sample": a.sample, "stage": "identify", "species": species,
               "confidence": conf, "method": method,
               "its_present": its_present, "its_identity": its_pid,
               "unite_hit": unite, "gather_top": gather}, open(a.out_json, "w"), indent=2)
    print(f"[id_classify] {a.sample}: {species} ({conf}, {method}, ITS%={its_pid})")


if __name__ == "__main__":
    main()
