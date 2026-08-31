#!/usr/bin/env python3
"""Stage 11 — summarise fungiSMASH output: count BGCs by product type (NRPS, PKS,
terpene, RiPP, hybrid…). MIBiG-known vs novel is resolved later by BiG-SCAPE across
all isolates in the comparative layer."""
from __future__ import annotations
import argparse, glob, json, os
from collections import Counter


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
    a = ap.parse_args()
    types = parse_antismash(a.as_dir)
    json.dump({"sample": a.sample, "stage": "bgc", "n_clusters": sum(types.values()),
               "by_type": dict(types.most_common()),
               "clusters": [{"type": t, "count": c} for t, c in types.most_common()]},
              open(a.out, "w"), indent=2)
    print(f"[bgc_summary] {a.sample}: {sum(types.values())} BGCs")


if __name__ == "__main__":
    main()
