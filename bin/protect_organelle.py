#!/usr/bin/env python3
"""Put organelle contigs back after haplotig purging (W4.1 follow-up).

  protect_organelle.py --assembly asm.fasta --purged purged.fa --tblastn mito_hits.tsv
                       [--min-genes 2] [--max-len 250000] [--max-gc 0.40]
                       --out restored.fasta --json protect.json

purge_dups judges contigs by read depth relative to the nuclear peak. A mitochondrial genome sits
far above that peak, so purge_dups classifies it as collapsed duplication and drops it. On a real
A. fumigatus assembly it removed nine of the ten organelle-shaped contigs, including the 30.7 kb,
25.5 % GC mitogenome, which emptied every mitochondrial column downstream. The behaviour appeared
only once purge_dups was fixed to run at all: before that it exited 127 and purged nothing.

So: classify the PRE-purge assembly with the same evidence stage 04 uses (>= --min-genes distinct
core mitochondrial proteins by tblastn, length and GC bounds), and append any contig meeting it
that purging dropped. Nothing else is restored — genuine haplotigs stay purged.

get_seqs renames what it keeps (`contig_7` -> `contig_7_1`), so presence is matched on the name
with a trailing `_<n>` removed, not on an exact string.
"""
from __future__ import annotations
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mito_extract import read_fasta, genes_per_contig, classify  # noqa: E402

SUFFIX = re.compile(r"_\d+$")


def base_names(names):
    """Names as purge_dups' get_seqs may have rewritten them, reduced back to the original."""
    out = set()
    for n in names:
        out.add(n)
        out.add(SUFFIX.sub("", n))
    return out


def restore(assembly, purged, tblastn, min_genes=2, max_len=250000, max_gc=0.40):
    pre, post = read_fasta(assembly), read_fasta(purged)
    genes = genes_per_contig(tblastn)
    mito, _why = classify(pre, genes, max_len=max_len, min_genes=min_genes, max_gc=max_gc)
    kept = base_names(post)
    dropped = [c for c in mito if c not in kept]
    return pre, post, sorted(mito), dropped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True, help="assembly as it was BEFORE purging")
    ap.add_argument("--purged", required=True, help="purge_dups output")
    ap.add_argument("--tblastn", required=True, help="tblastn of core mitochondrial proteins vs --assembly")
    ap.add_argument("--min-genes", type=int, default=2)
    ap.add_argument("--max-len", type=int, default=250000)
    ap.add_argument("--max-gc", type=float, default=0.40)
    ap.add_argument("--out", required=True)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    pre, post, mito, dropped = restore(a.assembly, a.purged, a.tblastn,
                                       a.min_genes, a.max_len, a.max_gc)
    with open(a.out, "w") as fh:
        for name, seq in post.items():
            fh.write(">%s\n%s\n" % (name, seq))
        for name in dropped:
            fh.write(">%s\n%s\n" % (name, pre[name]))
    doc = {"organelle_contigs": mito, "restored": dropped, "n_restored": len(dropped),
           "n_contigs_purged_output": len(post), "n_contigs_final": len(post) + len(dropped)}
    if a.json:
        json.dump(doc, open(a.json, "w"), indent=2)
    print("[protect_organelle] %d organelle contig(s) identified, %d restored after purging%s"
          % (len(mito), len(dropped), (": " + ", ".join(dropped)) if dropped else ""))


if __name__ == "__main__":
    main()
