#!/usr/bin/env python3
"""Read statistics for the read-QC stage JSON (W2.2): fastp's JSON for Illumina, NanoPlot's
NanoStats.txt for ONT (raw and filtered). Pure stdlib; every missing input yields nulls, never
a crash, and the stage JSON is written by merging into an existing document when present.

  read_stats.py --sample ID --platform illumina|ont|hybrid [--fastp fastp.json]
                [--nanostats-raw nanoplot_raw/rawNanoStats.txt] [--nanostats-filt ...]
                [--merge-into ID.readqc.json] --out ID.readqc.json
"""
from __future__ import annotations
import argparse, json, os, re


def fastp_stats(path):
    if not path or not os.path.exists(path):
        return None
    try:
        d = json.load(open(path))
    except Exception:  # noqa: BLE001
        return {"error": "fastp.json unreadable"}
    b, a = d.get("summary", {}).get("before_filtering", {}), d.get("summary", {}).get("after_filtering", {})
    return {"reads_raw": b.get("total_reads"), "bases_raw": b.get("total_bases"),
            "reads_trimmed": a.get("total_reads"), "bases_trimmed": a.get("total_bases"),
            "q30_rate_trimmed": a.get("q30_rate"), "gc_trimmed": a.get("gc_content"),
            "mean_len_trimmed": (a.get("total_bases") / a.get("total_reads")) if a.get("total_reads") else None,
            "duplication_rate": d.get("duplication", {}).get("rate"),
            "adapter_trimmed_reads": d.get("adapter_cutting", {}).get("adapter_trimmed_reads")}


_NANO_KEYS = {"Mean read length": "mean_len", "Mean read quality": "mean_q", "Median read length": "median_len",
              "Median read quality": "median_q", "Number of reads": "reads", "Read length N50": "n50",
              "Total bases": "bases"}


def nanostats(path):
    """NanoPlot NanoStats.txt -> dict; also handles the '>Q10' style rows."""
    if not path or not os.path.exists(path):
        return None
    out = {}
    for line in open(path):
        line = line.rstrip("\n")
        for key, name in _NANO_KEYS.items():
            if line.startswith(key):
                m = re.search(r"([-\d.,]+)\s*$", line)
                if m:
                    out[name] = float(m.group(1).replace(",", ""))
        m = re.match(r"^>Q(\d+):\s+([\d,]+)\s+\(([\d.]+)%\)", line)
        if m:
            out[f"reads_over_q{m.group(1)}_pct"] = float(m.group(3))
    for k in ("reads", "n50", "bases"):
        if k in out:
            out[k] = int(out[k])
    return out or None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--platform", required=True)
    ap.add_argument("--fastp"); ap.add_argument("--nanostats-raw"); ap.add_argument("--nanostats-filt")
    ap.add_argument("--merge-into"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    doc = {}
    if a.merge_into and os.path.exists(a.merge_into):
        try:
            doc = json.load(open(a.merge_into))
        except Exception:  # noqa: BLE001
            doc = {}
    doc.setdefault("sample", a.sample); doc.setdefault("stage", "readqc"); doc["platform"] = a.platform
    stats = {"illumina": fastp_stats(a.fastp), "ont_raw": nanostats(a.nanostats_raw), "ont_filtered": nanostats(a.nanostats_filt)}
    doc["read_stats"] = stats
    # headline numbers for the report
    il, on = stats["illumina"] or {}, stats["ont_filtered"] or {}
    doc["total_bases"] = (il.get("bases_trimmed") or 0) + (on.get("bases") or 0) or None
    doc["ont_n50"] = on.get("n50")
    json.dump(doc, open(a.out, "w"), indent=2)
    print(f"[read_stats] {a.sample} ({a.platform}): illumina {il.get('reads_trimmed', 'NA')} reads / {il.get('bases_trimmed', 'NA')} bp after trimming; "
          f"ONT {on.get('reads', 'NA')} reads / {on.get('bases', 'NA')} bp after filtering (N50 {on.get('n50', 'NA')})")


if __name__ == "__main__":
    main()
