#!/usr/bin/env python3
"""Stage 08b helper — species-aware BUSCO completeness -> JSON (W2.3).

  busco_lineage.py --sample ID --lineage eurotiales_odb10 --qc-lineage fungi_odb10
                   [--compleasm-dir compleasm --busco-dir busco | --skipped <reason>] --out ID.busco_lineage.json

Writes lineage, busco_complete (S+D from compleasm, or C from BUSCO), single/duplicated/
fragmented/missing when available, and `skipped` with the reason when no run was made."""
from __future__ import annotations
import argparse, glob, json, os, re

from assembly_qc import busco_complete  # noqa: E402  (bin/ is on sys.path / PATH)


def busco_breakdown(compleasm_dir, busco_dir):
    """S, D, F, M percentages from the compleasm summary or the BUSCO short summary."""
    for p in glob.glob(os.path.join(compleasm_dir or "", "summary.txt")):
        t = open(p).read()
        out = {}
        for key, tag in (("single", "S"), ("duplicated", "D"), ("fragmented", "F"), ("missing", "M")):
            m = re.search(rf"\b{tag}:\s*([0-9.]+)%", t)
            if m:
                out[key] = float(m.group(1))
        if out:
            return out
    for p in glob.glob(os.path.join(busco_dir or "", "short_summary*.txt")):
        t = open(p).read()
        m = re.search(r"C:\s*([0-9.]+)%\[S:\s*([0-9.]+)%,D:\s*([0-9.]+)%\],F:\s*([0-9.]+)%,M:\s*([0-9.]+)%", t)
        if m:
            return {"single": float(m.group(2)), "duplicated": float(m.group(3)), "fragmented": float(m.group(4)), "missing": float(m.group(5))}
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--lineage", required=True); ap.add_argument("--qc-lineage", default="fungi_odb10")
    ap.add_argument("--compleasm-dir", default=None); ap.add_argument("--busco-dir", default=None)
    ap.add_argument("--skipped", default=None, help="reason when BUSCO was not re-run")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    comp = None if a.skipped else busco_complete(a.compleasm_dir, a.busco_dir)
    doc = {"sample": a.sample, "stage": "busco_lineage", "lineage": a.lineage or None, "qc_lineage": a.qc_lineage,
           "busco_complete": comp, "skipped": a.skipped,
           **(busco_breakdown(a.compleasm_dir, a.busco_dir) if not a.skipped else {})}
    if not a.skipped and comp is None:
        doc["note"] = "BUSCO/compleasm produced no summary"
    json.dump(doc, open(a.out, "w"), indent=2)
    print(f"[busco_lineage] {a.sample}: {a.lineage} -> {comp if comp is not None else 'NA'}%" + (f" (skipped: {a.skipped})" if a.skipped else ""))


if __name__ == "__main__":
    main()
