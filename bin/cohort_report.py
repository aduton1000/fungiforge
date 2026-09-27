#!/usr/bin/env python3
"""Cohort summary: merged master table, HTML overview and MultiQC custom content (W3.3, stage 17).

  cohort_report.py --masters m1.tsv m2.tsv ... --outdir 04_summary [--cohort cohort.json]
                   [--provenance provenance.json] [--schema master_schema.json] [--version 0.2.0]

Writes <outdir>/master_fungi.tsv (one header, one row per isolate, schema order), validates it,
renders cohort_report.html (species and compartment breakdown, QC and gate outcomes, resistance
and mycotoxin tallies, assembly and annotation ranges, novelty, the cohort's clusters/tree/clonal
groups when stage 16 ran) and multiqc/fungiforge_mqc.json + multiqc/multiqc_config.yaml so
MultiQC renders the same tables when it is installed."""
from __future__ import annotations
import argparse, collections, csv, datetime, json, os, statistics, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_master import load_schema, validate  # noqa: E402

NA = "NA"


def read_masters(paths, columns):
    """Merge per-isolate master rows into schema order; later files win on duplicate samples."""
    rows, order = {}, []
    for p in paths:
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            continue
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                s = (r.get("sample") or "").strip()
                if not s:
                    continue
                if s not in rows:
                    order.append(s)
                rows[s] = r
    out = []
    for s in sorted(order):
        r = rows[s]
        out.append({c: (r.get(c) if r.get(c) not in (None, "") else NA) for c in columns})
    return out


def num(v):
    try:
        f = float(v)
        return f if f == f else None
    except (TypeError, ValueError):
        return None


def tally(rows, col, split=None, drop=("NA", "none", "")):
    c = collections.Counter()
    for r in rows:
        v = (r.get(col) or NA).strip()
        for part in (v.split(split) if split else [v]):
            part = part.strip()
            if part and part not in drop:
                c[part] += 1
    return c


def stats(rows, col):
    vals = [num(r.get(col)) for r in rows]
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    return {"n": len(vals), "min": min(vals), "median": statistics.median(vals), "max": max(vals),
            "mean": round(statistics.fmean(vals), 2)}


def summarise(rows, cohort=None, provenance=None):
    d = {"n_isolates": len(rows),
         "species": tally(rows, "species").most_common(),
         "compartments": tally(rows, "compartment").most_common(),
         "facilities": tally(rows, "facility").most_common(),
         "gates": tally(rows, "gate", drop=("",)).most_common(),
         "qc_pass": tally(rows, "qc_pass", drop=("",)).most_common(),
         "read_verdicts": tally(rows, "read_verdict", drop=("",)).most_common(),
         "novelty": tally(rows, "novelty", drop=("",)).most_common(),
         "resistant_classes": tally(rows, "resistant_classes", split=";").most_common(),
         "cyp51A_TR": tally(rows, "cyp51A_TR").most_common(),
         "mycotoxins": tally(rows, "mycotoxin_clusters", split=";").most_common(),
         "mating_types": tally(rows, "mating_type", drop=("", "NA")).most_common(),
         "stages_failed": tally(rows, "stages_failed", split=";").most_common(20),
         "id_flags": tally(rows, "id_flags", split=";").most_common(10),
         "metrics": {c: stats(rows, c) for c in ("assembly_len", "n_contigs", "busco_complete", "n_proteins", "n_bgc",
                                                 "te_percent", "coverage", "genome_size_est", "mito_size_kb", "n_cazymes")}}
    d["n_with_resistance"] = sum(1 for r in rows if (r.get("resistant_classes") or NA) not in (NA, "none", ""))
    d["n_gated"] = sum(1 for r in rows if (r.get("gate") or "").startswith("skipped"))
    d["n_stage_failures"] = sum(1 for r in rows if (r.get("stages_failed") or "none") not in ("none", NA, ""))
    if cohort:
        d["cohort"] = {k: cohort.get(k) for k in ("n_isolates", "species_clusters", "phylogenomics", "clusters_snp", "bgc_families") if k in cohort}
    if provenance:
        d["provenance"] = {"pipeline": (provenance.get("pipeline") or {}).get("version"),
                           "commit": (provenance.get("pipeline") or {}).get("commit_id"),
                           "run": (provenance.get("run") or {}).get("run_name"),
                           "started": (provenance.get("run") or {}).get("start")}
    return d


def multiqc_custom(summary, rows, outdir):
    """MultiQC custom content: a general-statistics table plus bar plots of the cohort tallies."""
    os.makedirs(outdir, exist_ok=True)
    gen = {}
    for r in rows:
        gen[r["sample"]] = {k: (num(r.get(k)) if num(r.get(k)) is not None else (r.get(k) or NA))
                            for k in ("species", "busco_complete", "assembly_len", "n_contigs", "n_proteins", "n_bgc",
                                      "te_percent", "coverage", "resistant_classes", "gate")}
    doc = {"id": "fungiforge", "section_name": "FungiForge cohort",
           "description": f"{summary['n_isolates']} isolates · {summary['n_with_resistance']} with a resistance call · {summary['n_gated']} stopped by a gate",
           "plot_type": "table", "pconfig": {"id": "fungiforge_master", "title": "FungiForge master table", "col1_header": "sample"},
           "data": gen}
    json.dump(doc, open(os.path.join(outdir, "fungiforge_mqc.json"), "w"), indent=2)
    for key, title in (("species", "Species"), ("resistant_classes", "Resistant classes"), ("mycotoxins", "Mycotoxin clusters"), ("novelty", "Novelty")):
        data = dict(summary[key])
        if not data:
            continue
        json.dump({"id": f"fungiforge_{key}", "section_name": f"FungiForge — {title}", "plot_type": "bargraph",
                   "pconfig": {"id": f"fungiforge_{key}_plot", "title": title, "ylab": "isolates"},
                   "data": {title: data}}, open(os.path.join(outdir, f"fungiforge_{key}_mqc.json"), "w"), indent=2)
    cfg = ["custom_logo_title: FungiForge", "title: FungiForge cohort report",
           f"subtitle: {summary['n_isolates']} isolates", "report_comment: >-",
           "    Generated by fungiforge; every value also lives in 04_summary/master_fungi.tsv and the per-isolate stage JSONs.",
           "show_analysis_paths: false", "custom_data:", "    fungiforge:", "        file_format: 'json'"]
    open(os.path.join(outdir, "multiqc_config.yaml"), "w").write("\n".join(cfg) + "\n")
    return os.path.join(outdir, "fungiforge_mqc.json")


def render(summary, rows, version, template_dir=None):
    ctx = {"s": summary, "rows": rows, "version": version,
           "generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}
    tdir = template_dir or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fungiforge", "report_templates")
    try:
        import jinja2
        if os.path.exists(os.path.join(tdir, "cohort_report.html.j2")):
            env = jinja2.Environment(loader=jinja2.FileSystemLoader(tdir), autoescape=True)
            return env.get_template("cohort_report.html.j2").render(**ctx)
    except Exception as e:  # noqa: BLE001
        sys.stderr.write(f"[cohort_report] template unavailable ({e}); writing the plain layout\n")
    import html as _h
    parts = [f"<!doctype html><meta charset=utf-8><title>FungiForge cohort</title><h1>FungiForge cohort — {summary['n_isolates']} isolates</h1>"]
    for key in ("species", "compartments", "gates", "resistant_classes", "mycotoxins", "novelty"):
        parts.append(f"<h2>{_h.escape(key)}</h2><ul>" + "".join(f"<li>{_h.escape(str(k))}: {v}</li>" for k, v in summary[key]) + "</ul>")
    return "".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--masters", nargs="+", required=True); ap.add_argument("--outdir", required=True)
    ap.add_argument("--cohort", default=None); ap.add_argument("--provenance", default=None)
    ap.add_argument("--schema", default=None); ap.add_argument("--version", default="0.2.0")
    ap.add_argument("--no-multiqc-export", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    schema = load_schema(a.schema)
    columns = [c["name"] for c in schema["columns"]]
    paths = []
    for m in a.masters:
        paths += sorted(os.path.join(m, f) for f in os.listdir(m) if f.endswith(".tsv")) if os.path.isdir(m) else [m]
    rows = read_masters(paths, columns)
    merged = os.path.join(a.outdir, "master_fungi.tsv")
    with open(merged, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    rep = validate(merged, schema)
    cohort = json.load(open(a.cohort)) if a.cohort and os.path.exists(a.cohort) else None
    prov = json.load(open(a.provenance)) if a.provenance and os.path.exists(a.provenance) else None
    summary = summarise(rows, cohort, prov)
    summary["schema_validation"] = {"valid": rep["valid"], "errors": rep["errors"][:20], "warnings": rep["warnings"], "schema_version": rep["schema_version"]}
    open(os.path.join(a.outdir, "cohort_report.html"), "w").write(render(summary, rows, a.version))
    json.dump(summary, open(os.path.join(a.outdir, "cohort_summary.json"), "w"), indent=2)
    if not a.no_multiqc_export:
        multiqc_custom(summary, rows, os.path.join(a.outdir, "multiqc"))
    print(f"[cohort_report] {len(rows)} isolates -> {merged} (schema {'valid' if rep['valid'] else 'INVALID: ' + str(len(rep['errors'])) + ' errors'}); "
          f"cohort_report.html, cohort_summary.json" + ("" if a.no_multiqc_export else ", multiqc/"))
    sys.exit(0 if rep["valid"] else 1)


if __name__ == "__main__":
    main()
