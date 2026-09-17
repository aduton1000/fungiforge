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
    # appended 0.2.0 (W2.1): "pass", or "skipped(<reason>)" when a gate stopped the isolate (read
    # triage: verdict non_fungal/human before assembly; assembly QC: verdict or qc_pass false)
    "gate",
    # appended 0.2.0 (W2.2): read-level triage verdict and the k-mer profile
    "read_verdict", "genome_size_est", "heterozygosity_pct", "ploidy_hint", "coverage",
    # appended 0.2.0 (W2.3): identification evidence — secondary loci that agree with ITS, flags
    # (tie/discordant/unavailable), MLST scheme:ST, species-aware BUSCO lineage and completeness
    "id_loci_agree", "id_flags", "mlst_st", "busco_lineage_specific", "busco_complete_specific",
    # appended 0.2.0 (W2.4): read-level support of the resistance calls and copy-number flags
    "resistance_read_support", "copy_number_flags",
    # appended 0.2.0 (W2.5): annotation size/quality and how the prediction was trained
    "n_proteins", "pct_pfam", "pct_go", "pct_eggnog", "pct_interpro", "annotation_training",
    # appended 0.2.0 (W2.6): extras counts (ploidy and mating_type above are now nQuire / Pfam-based)
    "n_secreted", "n_effectors", "n_cazymes", "n_phibase_hits",
    # appended 0.2.0 (W2.7): mitochondrial genome
    "mito_size_kb", "mito_core_genes", "mito_circular", "mito_copy_ratio", "mito_heteroplasmic_sites",
    # appended 0.2.0 (W2.8): mobile elements (te_percent and n_mycovirus above are now populated by stage 10)
    "te_ltr_pct", "n_genomad_virus", "n_genomad_plasmid", "n_mito_heg",
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


def read_support_summary(res):
    """confirmed:N;discordant:M;reads_only:K;insufficient:J, or assembly_only / NA."""
    rs = g(res, "summary", "read_support", default=None)
    if not isinstance(rs, dict):
        return "NA"
    if rs.get("mode") != "reads":
        return "assembly_only"
    return ";".join(f"{k}:{rs.get(k, 0)}" for k in ("confirmed", "discordant", "reads_only", "insufficient"))


def build_row(sample, compartment, facility, season, S):
    ident = S.get("identify", {})
    qc = S.get("assembly_qc", {})
    res = S.get("resistance", {})
    bgc = S.get("bgc", {})
    nov = S.get("novelty", {})
    ext = S.get("extras", {})
    ann = S.get("annotate", {})
    org = S.get("organelle", {})
    prd = S.get("predict", {})
    mob = S.get("mobile", {})
    pol = S.get("polish", {})
    dc  = S.get("decontam", {})
    tri = S.get("triage", {})
    km  = S.get("kmer", {})
    bl  = S.get("busco_lineage", {})
    ml  = g(ident, "mlst", default={})
    agree = g(ident, "concordance", "secondary_agree", default=[])
    flags = g(ident, "flags", default=[])
    tops = g(dc, "top_species", default=[])
    if not (isinstance(tops, list) and tops):
        tops = g(tri, "top_species", default=[])          # isolate stopped at read triage: use the read-level taxon
    top_taxon = (f"{tops[0].get('name', '?')} ({tops[0].get('pct', 'NA')}%)"
                 if isinstance(tops, list) and tops and isinstance(tops[0], dict) else "NA")
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
        "gate": f"skipped({S['gate'].get('reason', '?')})" if isinstance(S.get("gate"), dict) else "pass",
        "read_verdict": g(tri, "verdict", default="NA"),
        "genome_size_est": g(km, "genome_size_est", default="NA"),
        "heterozygosity_pct": g(km, "heterozygosity_pct", default="NA"),
        "ploidy_hint": g(km, "ploidy_hint", default="NA"),
        "coverage": g(km, "coverage_from_kmers", default=g(km, "coverage_from_read_bases", default="NA")),
        "id_loci_agree": ";".join(agree) if isinstance(agree, list) and agree else "none",
        "id_flags": ";".join(flags) if isinstance(flags, list) and flags else "none",
        "mlst_st": (f"{ml.get('scheme')}:ST{ml.get('st')}" if isinstance(ml, dict) and ml.get("scheme") and ml.get("st")
                    else f"{ml.get('scheme')}:ST-" if isinstance(ml, dict) and ml.get("scheme") else "NA"),
        "busco_lineage_specific": g(bl, "lineage", default="NA") if g(bl, "busco_complete", default=None) not in (None, "NA") else "NA",
        "busco_complete_specific": g(bl, "busco_complete", default="NA"),
        "resistance_read_support": read_support_summary(res),
        "copy_number_flags": ";".join(g(res, "summary", "copy_number_flags", default=[]) or []) or "none",
        "n_proteins": g(ann, "n_proteins"), "pct_pfam": g(ann, "pct_pfam"), "pct_go": g(ann, "pct_go"),
        "pct_eggnog": g(ann, "pct_eggnog"), "pct_interpro": g(ann, "pct_interpro"),
        "annotation_training": (f"{g(prd, 'augustus_species')}/{g(prd, 'busco_db')}/genemark:{g(prd, 'genemark')}"
                                if isinstance(prd, dict) and prd else "NA"),
        "n_secreted": g(ext, "n_secreted"), "n_effectors": g(ext, "n_effectors"),
        "n_cazymes": g(ext, "n_cazymes"), "n_phibase_hits": g(ext, "n_phibase_hits"),
        "mito_size_kb": round((org.get("mito_size_estimate") or org["mito_size"]) / 1000, 1) if isinstance(org, dict) and org.get("mito_size") else "NA",
        "mito_core_genes": f"{org['n_core_genes']}/15" if isinstance(org, dict) and "n_core_genes" in org and org.get("mito_present") else "NA",
        "mito_circular": g(org, "circular") if g(org, "mito_present", default=False) else "NA",
        "mito_copy_ratio": g(org, "copy_ratio"),
        "mito_heteroplasmic_sites": g(org, "heteroplasmy", "n_heteroplasmic_sites"),
        "te_ltr_pct": g(mob, "te_landscape", "ltr_pct"), "n_genomad_virus": g(mob, "n_genomad_virus"),
        "n_genomad_plasmid": g(mob, "n_genomad_plasmid"), "n_mito_heg": g(mob, "n_mito_heg"),
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
