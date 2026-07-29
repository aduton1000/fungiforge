#!/usr/bin/env python3
"""Summarise the RepeatModeler TE library: family counts by superfamily (LTR
Gypsy/Copia, TIR, LINE, Helitron, …) — the fungal 'mobile element' inventory."""
from __future__ import annotations
import argparse, json, re
from collections import Counter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--telib", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cls = Counter()
    try:
        for line in open(a.telib):
            if line.startswith(">"):
                m = re.search(r"#(\S+)", line)
                cls[m.group(1) if m else "Unknown"] += 1
    except Exception:
        pass
    total = sum(cls.values())
    json.dump({"stage": "te", "n_families": total,
               "by_class": dict(cls.most_common())}, open(a.out, "w"), indent=2)
    print(f"[te_summary] {total} TE families")


if __name__ == "__main__":
    main()
