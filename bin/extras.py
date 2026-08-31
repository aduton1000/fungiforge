#!/usr/bin/env python3
"""Stage 13 — eukaryote extras. Mating-type idiomorph (MAT1-1-1 / MAT1-2-1 by protein
name), a secretome/CAZyme/virulence tally when the annotation carries those terms, and
placeholders for ploidy/heterozygosity (nQuire/Smudgeplot run separately when reads are
available). Genome-only, degrades gracefully."""
from __future__ import annotations
import argparse, json, os, re

MAT = re.compile(r"MAT1-?1|MAT1-?2|alpha[- ]?box|HMG.?box.*mating", re.I)
CAZY = re.compile(r"glycoside hydrolase|glycosyltransferase|carbohydrate.?binding|GH\d+|GT\d+|AA\d+", re.I)
SECR = re.compile(r"signal peptide|secreted|effector", re.I)
VIR = re.compile(r"virulence|adhesin|melanin|phospholipase|candidalysin|siderophore", re.I)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--proteins", required=True)
    ap.add_argument("--ont"); ap.add_argument("--data-dir", default=""); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    mat, cazy, secr, vir = set(), 0, 0, 0
    try:
        for line in open(a.proteins):
            if line.startswith(">"):
                h = line
                if MAT.search(h):
                    m = MAT.search(h); mat.add(m.group(0))
                if CAZY.search(h): cazy += 1
                if SECR.search(h): secr += 1
                if VIR.search(h): vir += 1
    except Exception:
        pass
    mating = "MAT1-1" if any("1-1" in m or "alpha" in m.lower() for m in mat) else \
             ("MAT1-2" if any("1-2" in m or "hmg" in m.lower() for m in mat) else "undetermined")
    json.dump({"sample": a.sample, "stage": "extras", "mating_type": mating,
               "mating_hits": sorted(mat), "n_cazyme_like": cazy, "n_secreted_like": secr,
               "n_virulence_like": vir, "ploidy": "NA",
               "note": "ploidy/heterozygosity via nQuire/Smudgeplot when reads available; "
                       "secretome/CAZyme are annotation-keyword tallies (run_dbCAN/EffectorP refine)"},
              open(a.out, "w"), indent=2)
    print(f"[extras] {a.sample}: mating={mating}, cazy~{cazy}, vir~{vir}")


if __name__ == "__main__":
    main()
