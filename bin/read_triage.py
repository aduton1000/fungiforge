#!/usr/bin/env python3
"""Read-level triage stage JSON (W2.2, stage 00b): Kraken2 on a read subsample -> verdict.

  read_triage.py --sample ID --report k2.report --n-reads N --platform ont|illumina|hybrid
                 [--total-reads M] --out ID.triage.json

Reuses kraken_taxonomy.summarize logic (report-only lineage reconstruction) and the shared
verdict rules: non_fungal (bacteria+archaea+viruses >= 50 %), human (>= 50 %), fungal
(fungi >= 10 % and prokaryote < 20 %), likely_fungal (unclassified >= 50 % and prokaryote
< 10 %), else mixed. The same rules the contig-level verdict uses after assembly, so a
plate's bacteria are stopped before assembly (the gate) and the two verdicts are comparable.
"""
from __future__ import annotations
import argparse, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kraken_taxonomy as kt  # noqa: E402


def summarize_report(report):
    nodes = kt.parse_report(report)
    total = sum(n["taxon_n"] for n in nodes.values())
    groups = {}
    for t, n in nodes.items():
        g = kt.group_of(nodes, t)
        groups[g] = groups.get(g, 0) + n["taxon_n"]
    frac = {g: (c / total if total else 0.0) for g, c in groups.items()}
    return {"classified_total": total, "pct": {g: round(100 * v, 2) for g, v in frac.items()},
            "top_species": kt.top_species(nodes), "verdict": kt.verdict(frac) if total else "not_run"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--report", required=True)
    ap.add_argument("--n-reads", type=int, required=True, help="reads (or pairs) in the subsample classified")
    ap.add_argument("--total-reads", type=int, default=None); ap.add_argument("--platform", default="unknown")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = summarize_report(a.report) if os.path.exists(a.report) else {"classified_total": 0, "pct": {}, "top_species": [], "verdict": "not_run"}
    top = res["top_species"][0] if res["top_species"] else None
    doc = {"sample": a.sample, "stage": "triage", "platform": a.platform, "subsample_reads": a.n_reads,
           "total_reads": a.total_reads, "verdict": res["verdict"], "composition_pct": res["pct"],
           "top_species": res["top_species"],
           "top_taxon": f"{top['name']} ({top['pct']}%)" if top else None,
           "note": ("read-level Kraken2 verdict on a subsample; the PlusPF database holds few fungal genomes, "
                    "so a genuine fungus often reads as mostly unclassified with a fungal minority (likely_fungal)")}
    json.dump(doc, open(a.out, "w"), indent=2)
    print(f"[read_triage] {a.sample}: {doc['verdict']} on {a.n_reads} reads; composition {res['pct']}; top {doc['top_taxon']}")


if __name__ == "__main__":
    main()
