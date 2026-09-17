#!/usr/bin/env python3
"""Stage 11 — summarise fungiSMASH output: count BGCs by product type (NRPS, PKS, terpene,
RiPP, hybrid…) and, from antiSMASH's KnownClusterBlast results (W3.2, --cb-knownclusters),
record every region's best MIBiG hit and flag known mycotoxin / antibiotic clusters
(fungiforge/resources/mycotoxin_compounds.tsv). Cross-isolate families are built by
stage 16 (cohort) from the region GenBanks."""
from __future__ import annotations
import argparse, glob, json, os, re
from collections import Counter


def load_compounds(path):
    rows = []
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3 and f[0] != "name":
                    rows.append((f[0], f[1], re.compile(f[2], re.I)))
    return rows


def knowncluster_hits(as_dir):
    """Best MIBiG hit per region from the result JSON: {region_number: {accession, description,
    cluster_type, core_gene_hits, blast_score, n_hits, similarity}}; similarity = the fraction of
    the region's CDS with a hit in the MIBiG cluster (pairings / region CDS count)."""
    out = {}
    for jf in glob.glob(os.path.join(as_dir or "", "*.json")):
        try:
            data = json.load(open(jf))
        except Exception:  # noqa: BLE001
            continue
        for rec in data.get("records", []):
            n_cds = {}
            for feat in rec.get("features", []):
                if feat.get("type") == "region":
                    rn = int(feat.get("qualifiers", {}).get("region_number", ["0"])[0])
                    loc = feat.get("location", "")
                    n_cds[rn] = (loc, feat.get("qualifiers", {}).get("product", []))
            cb = rec.get("modules", {}).get("antismash.modules.clusterblast", {})
            kc = cb.get("knowncluster") or cb.get("knownclusterblast") or {}
            for res in kc.get("results", []):
                rn = res.get("region_number")
                rank = res.get("ranking") or []
                if not rank:
                    out[f"{rec.get('id', '?')}|{rn}"] = {"region_number": rn, "record": rec.get("id"), "n_hits": res.get("total_hits", 0), "best": None}
                    continue
                cluster, details = rank[0]
                pair = details.get("pairings") or []
                out[f"{rec.get('id', '?')}|{rn}"] = {
                    "region_number": rn, "record": rec.get("id"), "n_hits": res.get("total_hits", len(rank)),
                    "best": {"accession": cluster.get("accession"), "description": cluster.get("description"), "cluster_type": cluster.get("cluster_type"),
                             "core_gene_hits": details.get("core_gene_hits"), "blast_score": details.get("blast_score"),
                             "n_query_proteins_hit": len({p[0] for p in pair}) if pair else None, "n_reference_proteins": len(cluster.get("proteins") or [])},
                    "products": n_cds.get(rn, ("", []))[1]}
    return out


def flag_compounds(hits, compounds):
    flags = []
    for key, h in hits.items():
        b = h.get("best")
        if not b or not b.get("description"):
            continue
        for name, cls, rx in compounds:
            if rx.search(b["description"]):
                flags.append({"compound": name, "class": cls, "region": key, "mibig": b["accession"], "description": b["description"],
                              "core_gene_hits": b.get("core_gene_hits"), "n_query_proteins_hit": b.get("n_query_proteins_hit"), "n_reference_proteins": b.get("n_reference_proteins")})
                break
    return flags


def parse_antismash(as_dir):
    types = Counter()
    for jf in glob.glob(os.path.join(as_dir or "", "*.json")):
        try:
            data = json.load(open(jf))
        except Exception:
            continue
        for rec in data.get("records", []):
            for feat in rec.get("features", []):
                if feat.get("type") == "region":
                    for p in feat.get("qualifiers", {}).get("product", []):
                        types[p] += 1
    if not types:  # fallback: count region gbk files
        n = len(glob.glob(os.path.join(as_dir or "", "*region*.gbk")))
        if n: types["region"] = n
    return types


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--as-dir"); ap.add_argument("--out", required=True)
    ap.add_argument("--status", default=None, help="ok | failed | skipped (from the module)")
    ap.add_argument("--compounds", default=None, help="mycotoxin_compounds.tsv (default: packaged copy)")
    ap.add_argument("--log", default=None, help="antiSMASH log; its tail is recorded on failure")
    a = ap.parse_args()
    types = parse_antismash(a.as_dir)
    ran = bool(glob.glob(os.path.join(a.as_dir or "", "*.json")))   # antiSMASH wrote its result JSON
    status = a.status or ("ok" if ran else "failed")
    if status == "ok" and not ran:
        status = "failed"
    comp_path = a.compounds or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fungiforge", "resources", "mycotoxin_compounds.tsv")
    if not os.path.exists(comp_path):
        try:
            from importlib.resources import files as _pkg_files
            comp_path = str(_pkg_files("fungiforge").joinpath("resources", "mycotoxin_compounds.tsv"))
        except Exception:  # noqa: BLE001
            pass
    hits = knowncluster_hits(a.as_dir) if status == "ok" else {}
    flags = flag_compounds(hits, load_compounds(comp_path))
    known = [h for h in hits.values() if h.get("best")]
    out = {"sample": a.sample, "stage": "bgc", "status": status,
           "n_clusters": sum(types.values()) if status == "ok" else None,
           "by_type": dict(types.most_common()),
           "clusters": [{"type": t, "count": c} for t, c in types.most_common()],
           "knownclusterblast": {"ran": bool(hits), "n_regions_with_mibig_hit": len(known),
                                 "best_hits": {k: v["best"] for k, v in hits.items() if v.get("best")}},
           "mycotoxin_flags": flags,
           "mycotoxin_compounds": sorted({f["compound"] for f in flags if f["class"] == "mycotoxin"}),
           "bioactive_compounds": sorted({f["compound"] for f in flags if f["class"] != "mycotoxin"})}
    if status == "ok" and not hits:
        out["knownclusterblast"]["note"] = "no KnownClusterBlast results in the antiSMASH JSON (run with --cb-knownclusters and the knownclusterblast database)"
    if status != "ok":
        tail = ""
        if a.log and os.path.exists(a.log):
            tail = "".join(open(a.log, errors="ignore").readlines()[-8:]).strip()
        out["note"] = f"antiSMASH {status}" + (f": {tail}" if tail else "")
    json.dump(out, open(a.out, "w"), indent=2)
    print(f"[bgc_summary] {a.sample}: {status}; {out['n_clusters']} BGCs; MIBiG hits {out['knownclusterblast']['n_regions_with_mibig_hit']}; "
          f"mycotoxins {','.join(out['mycotoxin_compounds']) or 'none'}")


if __name__ == "__main__":
    main()
