#!/usr/bin/env python3
"""Mycovirus / endogenous viral element screen from DIAMOND blastx hits against RVDB-prot (W2.8).

  mycovirus_screen.py --sample ID --hits rvdb.tsv --out mycovirus.json
                      [--min-pident 35] [--min-alen 80] [--mito-contigs c1,c2]

rvdb.tsv: diamond blastx outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle"
of the assembly contigs (queries) against RVDB-prot. Hits are grouped into loci per contig and
classified by the virus name in the RVDB title: MYCOVIRUS families (dsRNA/ssRNA fungal viruses),
RETROELEMENT-like (LTR retrotransposons, retrovirus-like: RepeatModeler's business, not a virus),
and OTHER. WGS is DNA-only: RNA mycoviruses are only visible when reverse-transcribed into the
genome (endogenous viral elements) — stated in the output."""
from __future__ import annotations
import argparse, collections, json, os, re

MYCO = re.compile(r"partitivir|totivir|victorivir|chrysovir|hypovir|narnavir|mitovir|endornavir|megabirnavir|quadrivir|botrexvir|alphaflexivir|"
                  r"polymycovir|fusarivir|yadokarivir|ourmiavir|botourmiavir|alternavir|umbra-?like|mycovir|mycoflexivir|deltaflexivir|"
                  r"gammaflexivir|hadakavir|curvularivir|fusagravir|phlegivir|amalgavir|bunyavir.*(fung|mycol)|tymo-?like|mycotymovir|"
                  r"mycoreovir|reovir.*(fung|rosellinia)|fungal (virus|dsRNA)|(fusarium|aspergillus|penicillium|botrytis|sclerotinia|rhizoctonia|cryphonectria|"
                  r"magnaporthe|colletotrichum|alternaria|trichoderma|rosellinia|hypoxylon|phomopsis|diaporthe|leptosphaeria|talaromyces|beauveria|"
                  r"metarhizium|ustilago|agaricus|lentinula|pleurotus|cryptococcus|saccharomyces|candida) [a-z ]*virus", re.I)
RETRO = re.compile(r"retrotranspos|retrovir|\bty[13]\b|gypsy|copia|metavir|pseudovir|belpaovir|caulimovir|\bLTR\b|reverse transcriptase.*(retro|element)|"
                   r"endogenous retrovirus|hepadnavir|transposon|polinton|maverick", re.I)


def parse_hits(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            try:
                qs, qe = int(f[4]), int(f[5])
                rows.append({"contig": f[0], "subject": f[1], "pident": float(f[2]), "alen": int(f[3]), "qstart": min(qs, qe), "qend": max(qs, qe),
                             "strand": "+" if qe >= qs else "-", "evalue": float(f[8]), "bitscore": float(f[9]), "title": f[10]})
            except ValueError:
                continue
    return rows


def classify_title(title):
    t = title.split("|")[-1] if "|" in title else title
    if RETRO.search(t):
        return "retroelement_like"
    if MYCO.search(t):
        return "mycovirus"
    return "other_viral"


def virus_name(title):
    m = re.search(r"\[([^\]]+)\]\s*$", title)
    return m.group(1) if m else title.split("|")[-1][:80]


def loci(rows, min_pident=35.0, min_alen=80, merge_gap=1000):
    """Merge overlapping/nearby hits per contig into loci with the best hit's class and name."""
    keep = sorted((r for r in rows if r["pident"] >= min_pident and r["alen"] >= min_alen), key=lambda r: (r["contig"], r["qstart"]))
    out = []
    for r in keep:
        if out and out[-1]["contig"] == r["contig"] and r["qstart"] <= out[-1]["end"] + merge_gap:
            L = out[-1]
            L["end"] = max(L["end"], r["qend"]); L["n_hits"] += 1
            if r["bitscore"] > L["bitscore"]:
                L.update(bitscore=r["bitscore"], pident=r["pident"], best=virus_name(r["title"]), category=classify_title(r["title"]), subject=r["subject"])
        else:
            out.append({"contig": r["contig"], "start": r["qstart"], "end": r["qend"], "strand": r["strand"], "n_hits": 1, "bitscore": r["bitscore"],
                        "pident": r["pident"], "best": virus_name(r["title"]), "category": classify_title(r["title"]), "subject": r["subject"]})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--hits", default=None); ap.add_argument("--out", required=True)
    ap.add_argument("--min-pident", type=float, default=35.0); ap.add_argument("--min-alen", type=int, default=80)
    ap.add_argument("--mito-contigs", default="")
    a = ap.parse_args()
    mito = {c for c in a.mito_contigs.split(",") if c}
    rows = parse_hits(a.hits)
    L = loci(rows, a.min_pident, a.min_alen)
    for x in L:
        x["compartment"] = "mitochondrial" if x["contig"] in mito else "nuclear"
    by = collections.Counter(x["category"] for x in L)
    myco = [x for x in L if x["category"] == "mycovirus"]
    doc = {"sample": a.sample, "n_hits": len(rows), "n_loci": len(L), "by_category": dict(by),
           "n_mycovirus_eve_candidates": len(myco), "mycovirus_candidates": sorted(myco, key=lambda x: -x["bitscore"])[:100],
           "retroelement_like_loci": by.get("retroelement_like", 0), "other_viral_loci": by.get("other_viral", 0),
           "families": dict(collections.Counter(x["best"] for x in myco).most_common(30)),
           "note": "DNA-only: RNA mycoviruses are visible only as endogenous viral elements (reverse-transcribed copies); "
                   "retroelement-like hits are LTR retrotransposons, reported separately; confirm candidates with RNA-seq / dsRNA extraction"}
    json.dump(doc, open(a.out, "w"), indent=2)
    print(f"[mycovirus_screen] {a.sample}: {len(rows)} hits -> {len(L)} loci: mycovirus-like {len(myco)}, retroelement-like {doc['retroelement_like_loci']}, other {doc['other_viral_loci']}")


if __name__ == "__main__":
    main()
