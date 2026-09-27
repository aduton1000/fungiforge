#!/usr/bin/env python3
"""Derive resistance-panel rows from the FungAMR catalogue (W2.4).

FungAMR (Landry lab, Nat Microbiol 2025) lists every published resistance mutation with a
confidence score (1 = engineered in the same species and measured … 7 = population
association, 8 = seen in a resistant natural isolate without validation; negative = shown not
to cause resistance). This module turns the table into panel rows compatible with
fungiforge/resources/af_resistance_panel.tsv, one per (species, gene): simple substitutions
(e.g. Y132F) become known_mutations with their strongest positive tier; TR/deletion/
overexpression tokens are left to the curated overlay. Species names are normalised through
the catalogue's own synonyms so a call of "Candida auris" matches "Candidozyma auris".

  fungamr_panel.py --data-dir <db> [--out <db>/fungamr/fungamr_panel.tsv] [--min-tier 8]
"""
from __future__ import annotations
import argparse, csv, os, re, sys
from collections import defaultdict

SUB_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")
DEL_RE = re.compile(r"^([A-Z])(\d+)(del|\*|X)$", re.I)
SYNONYMS = {  # FungAMR name -> the names the identification stage may produce
    "Candidozyma auris": ["Candida auris"],
    "Nakaseomyces glabratus": ["Candida glabrata", "Nakaseomyces glabrata"],
    "Pichia kudriavzevii": ["Candida krusei"],
    "Candidozyma haemuli": ["Candida haemulonii"],
    "Clavispora lusitaniae": ["Candida lusitaniae"],
    "Meyerozyma guilliermondii": ["Candida guilliermondii"],
    "Kluyveromyces marxianus": ["Candida kefyr"],
    "Nakaseomyces bracarensis": ["Candida bracarensis"],
    "Nakaseomyces nivariensis": ["Candida nivariensis"],
}
DEFAULT_CLASS = {"fluconazole": "azole", "voriconazole": "azole", "itraconazole": "azole", "posaconazole": "azole",
                 "isavuconazole": "azole", "clotrimazole": "azole", "miconazole": "azole", "ketoconazole": "azole",
                 "caspofungin": "echinocandin", "micafungin": "echinocandin", "anidulafungin": "echinocandin",
                 "rezafungin": "echinocandin", "ibrexafungerp": "triterpenoid", "amphotericin b": "polyene",
                 "nystatin": "polyene", "5-fluorocytosine": "pyrimidine", "flucytosine": "pyrimidine",
                 "terbinafine": "allylamine", "olorofim": "orotomide", "fosmanogepix": "gpi_inhibitor",
                 "manogepix": "gpi_inhibitor"}


def drug_classes(path):
    """drugs_class.csv (drug;usage;class) -> {drug_lower: class_slug}; falls back to DEFAULT_CLASS."""
    out = dict(DEFAULT_CLASS)
    if path and os.path.exists(path):
        with open(path, encoding="utf-8", errors="replace") as fh:
            rd = csv.reader(fh, delimiter=";")
            next(rd, None)
            for row in rd:
                if len(row) >= 3 and row[0].strip():
                    cls = row[2].strip().lower()
                    slug = {"azoles": "azole", "clinical azoles": "azole", "agricultural azoles": "azole_agricultural",
                            "echinocandins": "echinocandin", "polyenes": "polyene", "pyrimidine analogs": "pyrimidine",
                            "pyrimidine analogues": "pyrimidine", "allylamines": "allylamine"}.get(cls, re.sub(r"[^a-z0-9]+", "_", cls).strip("_"))
                    out.setdefault(row[0].strip().lower(), slug)
    return out


def load_rows(data_dir):
    tsv = None
    for cand in sorted(os.listdir(os.path.join(data_dir, "fungamr"))) if os.path.isdir(os.path.join(data_dir, "fungamr")) else []:
        if re.match(r"FungAMR_\d+\.tsv$", cand):
            tsv = os.path.join(data_dir, "fungamr", cand)
    if not tsv:
        return None, None
    with open(tsv, encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows, os.path.basename(tsv)


def derive_panel(rows, classes, min_tier=8):
    """-> list of panel rows (dict with the curated columns + evidence per mutation) and the
    evidence map {(gene_lower, species_lower, mutation): tier}."""
    per = defaultdict(lambda: {"muts": {}, "mut_classes": defaultdict(set), "drugs": set(), "classes": set(), "pmids": set()})
    for r in rows:
        try:
            tier = float((r.get("confidence score") or "").strip())
        except ValueError:
            continue
        if tier != tier or tier <= 0 or tier > min_tier:   # NaN, negative (shown not to confer) or beyond the tier cut
            continue
        gene = (r.get("gene or protein name") or "").strip()
        sp = (r.get("species") or "").strip()
        if not gene or not sp:
            continue
        drug = (r.get("drug") or "").strip()
        cls = classes.get(drug.lower())
        for tok in re.split(r"[,;/]+", (r.get("mutation") or "").strip()):
            tok = tok.strip()
            m = SUB_RE.match(tok)
            if not m or m.group(1) == m.group(3):
                continue
            key = (gene, sp)
            prev = per[key]["muts"].get(tok)
            per[key]["muts"][tok] = int(tier) if prev is None else min(prev, int(tier))
            if cls:
                per[key]["mut_classes"][tok].add(cls)
            if drug:
                per[key]["drugs"].add(drug.lower())
            if cls:
                per[key]["classes"].add(cls)
            pm = (r.get("pubmedid") or "").strip()
            if pm:
                per[key]["pmids"].add(pm)
    panel, evidence = [], {}
    for (gene, sp), d in sorted(per.items()):
        if not d["muts"]:
            continue
        names = [sp] + SYNONYMS.get(sp, [])
        muts = sorted(d["muts"], key=lambda t: (int(SUB_RE.match(t).group(2)), t))
        positions = sorted({int(SUB_RE.match(t).group(2)) for t in muts})
        row = {"gene": gene, "organism_regex": "|".join(re.escape(n) for n in names), "organism": sp,
               "drug_class": ";".join(sorted(d["classes"])) or "unknown", "drugs": ",".join(sorted(d["drugs"])),
               "mechanism": "substitution", "hotspot_aa": ",".join(str(p) for p in positions),
               "known_mutations": ",".join(muts), "tiers": ",".join(f"{t}:{d['muts'][t]}" for t in muts),
               "mut_classes": ",".join(f"{t}:{'/'.join(sorted(d['mut_classes'][t]))}" for t in muts),
               "note": f"FungAMR-derived; {len(muts)} substitutions, best tier {min(d['muts'].values())}",
               "source": "FungAMR " + ";".join(sorted(d["pmids"])[:5])}
        panel.append(row)
        for t, tier in d["muts"].items():
            for n in names:
                evidence[(gene.lower(), n.lower(), t)] = {"tier": tier, "classes": sorted(d["mut_classes"][t])}
    return panel, evidence


def write_panel(panel, path):
    cols = ["gene", "organism_regex", "drug_class", "drugs", "mechanism", "hotspot_aa", "known_mutations", "tiers", "mut_classes", "note", "source"]
    with open(path, "w") as fh:
        fh.write("# FungAMR-derived resistance panel (bin/fungamr_panel.py). Columns as af_resistance_panel.tsv plus\n"
                 "# `tiers` (mutation:FungAMR confidence score, 1 strongest .. 8 unvalidated natural isolate) and\n# `mut_classes` (mutation:drug classes the score was recorded for).\n")
        fh.write("\t".join(cols) + "\n")
        for r in panel:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def load_or_derive(data_dir, min_tier=8):
    """Panel rows + evidence map from <data_dir>/fungamr: the pre-built fungamr_panel.tsv when
    present, else derived from the raw table on the fly (sub-second). ([], {}) without FungAMR."""
    if not data_dir:
        return [], {}, None
    fdir = os.path.join(data_dir, "fungamr")
    pre = os.path.join(fdir, "fungamr_panel.tsv")
    if os.path.exists(pre):
        panel, evidence = [], {}
        with open(pre) as fh:
            header = None
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    header = parts; continue
                row = dict(zip(header, parts))
                panel.append(row)
                names = [n for n in re.split(r"\|", row.get("organism_regex", "")) if n]
                mcls = {}
                for kv in (row.get("mut_classes") or "").split(","):
                    if ":" in kv:
                        t, c = kv.split(":", 1)
                        mcls[t] = [x for x in c.split("/") if x]
                for kv in (row.get("tiers") or "").split(","):
                    if ":" in kv:
                        t, tier = kv.split(":", 1)
                        for n in names:
                            evidence[(row["gene"].lower(), re.sub(r"\\", "", n).lower(), t)] = {"tier": int(tier), "classes": mcls.get(t, [])}
        return panel, evidence, "fungamr_panel.tsv"
    rows, src = load_rows(data_dir)
    if not rows:
        return [], {}, None
    panel, evidence = derive_panel(rows, drug_classes(os.path.join(fdir, "drugs_class.csv")), min_tier)
    return panel, evidence, src


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True); ap.add_argument("--out", default=None); ap.add_argument("--min-tier", type=int, default=8)
    a = ap.parse_args()
    rows, src = load_rows(a.data_dir)
    if not rows:
        sys.exit(f"no FungAMR_*.tsv under {a.data_dir}/fungamr")
    panel, evidence = derive_panel(rows, drug_classes(os.path.join(a.data_dir, "fungamr", "drugs_class.csv")), a.min_tier)
    out = a.out or os.path.join(a.data_dir, "fungamr", "fungamr_panel.tsv")
    write_panel(panel, out)
    print(f"[fungamr_panel] {src}: {len(panel)} species/gene rows, {len(evidence)} mutation entries -> {out}")


if __name__ == "__main__":
    main()
