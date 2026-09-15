#!/usr/bin/env python3
"""Stage 14 — merge every per-stage result.json into (a) a self-contained per-isolate
HTML report and (b) one row of the master table (the Layer-2 handoff contract).

The master TSV schema is fixed and organism-agnostic so the comparative analysis
layer (analysis/) can bind isolates to air/human/surface metadata exactly like the
bacterial `master_results` table. Missing stages -> NA (never a hard failure).
"""
from __future__ import annotations
import argparse
import glob
import html
import json
import os

# Fixed master-table columns (the Layer-2 contract). Keep additions append-only.
MASTER_COLS = [
    "sample", "compartment", "facility", "season",
    "species", "species_confidence", "id_method",
    "qc_pass", "busco_complete", "busco_lineage", "assembly_len", "n_contigs",
    "ploidy", "mating_type",
    "resistant_classes", "n_known_af_mutations", "cyp51A_TR",
    "novelty", "n_bgc", "n_mycovirus", "te_percent",
    "polish_mode",
    # appended 0.1.1: Stage 04 Kraken2 verdict (fungal | likely_fungal | non_fungal | human | mixed | not_run)
    "sample_verdict", "contam_removed_pct", "top_taxon",
    # appended 0.2.0: stage status contract — any stage not "ok" as stage:status (or "none")
    "stages_failed",
]


def load_jsons(paths):
    by_stage = {}
    for p in paths:
        try:
            d = json.load(open(p))
        except Exception:
            continue
        by_stage[d.get("stage", os.path.basename(p))] = d
    return by_stage


def g(d, *keys, default="NA"):
    """Nested get with default."""
    cur = d
    for k in keys:
        if isinstance(cur, dict) and k in cur:
            cur = cur[k]
        else:
            return default
    return cur if cur not in (None, "", []) else default


def build_row(sample, compartment, facility, season, S):
    ident = S.get("identify", {})
    qc = S.get("assembly_qc", {})
    res = S.get("resistance", {})
    bgc = S.get("bgc", {})
    nov = S.get("novelty", {})
    ext = S.get("extras", {})
    mob = S.get("mobile", {})
    pol = S.get("polish", {})
    dc  = S.get("decontam", {})
    tops = g(dc, "top_species", default=[])
    top_taxon = f"{tops[0]['name']} ({tops[0]['pct']}%)" if isinstance(tops, list) and tops else "NA"
    rc = g(res, "summary", "resistant_drug_classes", default=[])
    tr = g(res, "cyp51A_TR", "tr_type", default="NA")
    row = {
        "sample": sample, "compartment": compartment, "facility": facility, "season": season,
        "species": g(ident, "species"), "species_confidence": g(ident, "confidence"),
        "id_method": g(ident, "method"),
        "qc_pass": g(qc, "qc_pass"), "busco_complete": g(qc, "busco_complete"),
        "busco_lineage": g(qc, "lineage"), "assembly_len": g(qc, "assembly_len"),
        "n_contigs": g(qc, "n_contigs"),
        "ploidy": g(ext, "ploidy"), "mating_type": g(ext, "mating_type"),
        "resistant_classes": ";".join(rc) if isinstance(rc, list) and rc else "none",
        "n_known_af_mutations": g(res, "summary", "n_known_mutations", default=0),
        "cyp51A_TR": tr if tr else "NA",
        "novelty": g(nov, "verdict", default=g(nov, "novelty", default="NA")),
        "n_bgc": g(bgc, "n_clusters", default="NA") if g(bgc, "status", default="ok") == "ok" else "NA",
        "n_mycovirus": g(mob, "n_mycovirus", default="NA"),
        "te_percent": g(mob, "te_percent", default="NA"),
        "polish_mode": g(pol, "mode"),
        "sample_verdict": g(dc, "verdict", default="NA"),
        "contam_removed_pct": g(dc, "dropped_bp_pct", default="NA"),
        "top_taxon": top_taxon,
        "stages_failed": ";".join(f"{st}:{d.get('status')}" for st, d in sorted(S.items())
                                  if d.get("status") not in (None, "ok")) or "none",
    }
    return row


def render_html(sample, meta, S, row):
    def esc(x): return html.escape(str(x))
    rows = "".join(f"<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>" for k, v in row.items())
    # resistance calls table
    calls = S.get("resistance", {}).get("calls", [])
    call_rows = "".join(
        f"<tr><td>{esc(c.get('gene'))}</td><td>{esc(c.get('drug_class'))}</td>"
        f"<td>{esc(c.get('change', c.get('status')))}</td><td>{esc(c.get('known'))}</td>"
        f"<td>{esc(c.get('confidence'))}</td></tr>" for c in calls) or "<tr><td colspan=5>none</td></tr>"
    stages = ", ".join(sorted(S.keys()))
    try:
        import jinja2  # noqa: F401 — richer templating available if desired
    except Exception:
        pass
    return f"""<!doctype html><meta charset=utf-8><title>FungiForge · {esc(sample)}</title>
<style>body{{font:14px/1.5 system-ui,sans-serif;max-width:900px;margin:2rem auto;padding:0 1rem;color:#1a2a44}}
h1{{color:#17233d}}table{{border-collapse:collapse;width:100%;margin:1rem 0}}
th,td{{border:1px solid #dcdfe6;padding:6px 10px;text-align:left}}th{{background:#f4f6fa;width:34%}}
.k{{color:#8a6d1f}}caption{{text-align:left;font-weight:bold;margin:.5rem 0}}</style>
<h1>FungiForge — {esc(sample)}</h1>
<p><b>{esc(row['species'])}</b> · {esc(meta['compartment'])} / {esc(meta['facility'])} / {esc(meta['season'])}
· stages: {esc(stages)}</p>
<table><caption>Summary (master row)</caption>{rows}</table>
<table><caption>Antifungal-resistance calls</caption>
<tr><th>gene</th><th>class</th><th>change/status</th><th>known</th><th>confidence</th></tr>{call_rows}</table>
<p style="color:#777">Resistance calls on ONT-only assemblies are <b>provisional</b> until hybrid-polished; hybrid and Illumina-only assemblies are high-confidence (no homopolymer-indel risk).</p>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--compartment", default="NA"); ap.add_argument("--facility", default="NA"); ap.add_argument("--season", default="NA")
    ap.add_argument("--jsons", nargs="+", required=True)
    ap.add_argument("--html", required=True); ap.add_argument("--master", required=True)
    a = ap.parse_args()
    # expand any globs / dirs
    paths = []
    for j in a.jsons:
        paths += glob.glob(j) if any(c in j for c in "*?[") else [j]
    S = load_jsons(paths)
    meta = {"compartment": a.compartment, "facility": a.facility, "season": a.season}
    row = build_row(a.sample, a.compartment, a.facility, a.season, S)
    with open(a.master, "w") as fh:
        fh.write("\t".join(MASTER_COLS) + "\n")
        fh.write("\t".join(str(row.get(c, "NA")) for c in MASTER_COLS) + "\n")
    open(a.html, "w").write(render_html(a.sample, meta, S, row))
    print(f"[make_report] {a.sample}: {row['species']} · {len(S)} stages · report+master written")


if __name__ == "__main__":
    main()
