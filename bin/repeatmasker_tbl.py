#!/usr/bin/env python3
"""Parse a RepeatMasker `.tbl` summary into te_percent and a class breakdown (W2.8, stage 10).

  repeatmasker_tbl.py --tbl rm/<asm>.tbl --out te_landscape.json

Reads the header (sequences, total length, GC, bases masked) and every class/superfamily row
("LTR elements:  712  877221 bp  2.93 %"), keeping the indentation as the hierarchy level."""
from __future__ import annotations
import argparse, json, os, re

ROW = re.compile(r"^(\s*)([A-Za-z0-9/\-\.'() ,]+?):?\s+(\d+)\s+(\d+) bp\s+([0-9.]+) %\s*$")
ROW_NOCOUNT = re.compile(r"^(\s*)([A-Za-z0-9/\-\.'() ,]+?):\s+(\d+) bp\s+([0-9.]+) %\s*$")   # "Total interspersed repeats:  1207448 bp  4.03 %"


def parse_tbl(path):
    if not path or not os.path.exists(path):
        return None
    doc = {"classes": {}, "rows": []}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^sequences:\s+(\d+)", line)
            if m:
                doc["n_sequences"] = int(m.group(1)); continue
            m = re.match(r"^total length:\s+(\d+) bp", line)
            if m:
                doc["total_length"] = int(m.group(1)); continue
            m = re.match(r"^GC level:\s+([0-9.]+) %", line)
            if m:
                doc["gc_pct"] = float(m.group(1)); continue
            m = re.match(r"^bases masked:\s+(\d+) bp \(\s*([0-9.]+) %\)", line)
            if m:
                doc["bases_masked"] = int(m.group(1)); doc["te_percent"] = float(m.group(2)); continue
            m = ROW.match(line.rstrip("\n"))
            m2 = None if m else ROW_NOCOUNT.match(line.rstrip("\n"))
            if m or m2:
                if m:
                    indent, name, n, bp, pct = len(m.group(1)), m.group(2).strip(), int(m.group(3)), int(m.group(4)), float(m.group(5))
                else:
                    indent, name, n, bp, pct = len(m2.group(1)), m2.group(2).strip(), None, int(m2.group(3)), float(m2.group(4))
                doc["rows"].append({"name": name, "level": indent // 2, "n": n, "bp": bp, "pct": pct})
                if indent == 0 or name in ("SINEs", "LINEs", "LTR elements", "Ty1/Copia", "Gypsy/DIRS1", "Penelope", "hobo-Activator", "Tc1-IS630-Pogo",
                                            "En-Spm", "MULE-MuDR", "PiggyBac", "Tourist/Harbinger", "Rolling-circles", "Unclassified", "Total interspersed repeats",
                                            "Small RNA", "Satellites", "Simple repeats", "Low complexity"):
                    doc["classes"][name] = {"n": n, "bp": bp, "pct": pct}
    if "te_percent" not in doc and "Total interspersed repeats" in doc["classes"]:
        doc["te_percent"] = doc["classes"]["Total interspersed repeats"]["pct"]
    c = doc["classes"]
    doc["summary"] = {"te_percent": doc.get("te_percent"), "interspersed_pct": c.get("Total interspersed repeats", {}).get("pct"),
                      "ltr_pct": c.get("LTR elements", {}).get("pct"), "gypsy_pct": c.get("Gypsy/DIRS1", {}).get("pct"), "copia_pct": c.get("Ty1/Copia", {}).get("pct"),
                      "line_pct": c.get("LINEs", {}).get("pct"), "dna_transposon_pct": c.get("DNA transposons", {}).get("pct"),
                      "helitron_pct": c.get("Rolling-circles", {}).get("pct"), "unclassified_pct": c.get("Unclassified", {}).get("pct"),
                      "simple_repeat_pct": c.get("Simple repeats", {}).get("pct")}
    return doc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tbl", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = parse_tbl(a.tbl)
    json.dump(d or {"te_percent": None, "note": "no RepeatMasker table"}, open(a.out, "w"), indent=2)
    print(f"[repeatmasker_tbl] te_percent={d.get('te_percent') if d else 'NA'} ({len(d['rows']) if d else 0} rows)")


if __name__ == "__main__":
    main()
