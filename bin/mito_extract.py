#!/usr/bin/env python3
"""Separate the mitochondrial genome from an assembly (W2.7, stage 04).

  mito_extract.py --sample ID --assembly asm.fa --tblastn hits.tsv --out-mito ID.mito.fasta
                  --out-nuclear ID.nuclear.fasta --json ID.mito_extract.json [--max-len 250000]
                  [--min-genes 2] [--max-gc 0.40]

hits.tsv: tblastn of fungiforge/resources/mito/mito_core_proteins.faa against the assembly,
outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore". A contig is
mitochondrial when it carries at least --min-genes distinct core genes (E <= 1e-10), is at most
--max-len long and has a GC fraction below --max-gc (fungal mtDNA is AT-rich; the GC guard keeps a
nuclear NUMT-carrying chromosome out). Everything else stays nuclear."""
from __future__ import annotations
import argparse, json, os


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def gc(seq):
    s = seq.upper()
    n = sum(1 for c in s if c in "ACGT")
    return round(sum(1 for c in s if c in "GC") / n, 3) if n else 0.0


def genes_per_contig(tblastn_path, max_evalue=1e-10):
    per = {}
    if not tblastn_path or not os.path.exists(tblastn_path):
        return per
    with open(tblastn_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            try:
                ev = float(f[8])
            except ValueError:
                continue
            if ev <= max_evalue:
                per.setdefault(f[1], set()).add(f[0].split("|")[0])
    return per


def classify(seqs, genes, max_len=250000, min_genes=2, max_gc=0.40):
    mito, why = [], {}
    for name, seq in seqs.items():
        g = sorted(genes.get(name, set()))
        gcf = gc(seq)
        if len(g) >= min_genes and len(seq) <= max_len and gcf < max_gc:
            mito.append(name)
            why[name] = {"length": len(seq), "gc": gcf, "core_genes": g, "mitochondrial": True}
        elif g:
            why[name] = {"length": len(seq), "gc": gcf, "core_genes": g, "mitochondrial": False,
                         "reason": "too long" if len(seq) > max_len else "GC too high (NUMT?)" if gcf >= max_gc else "fewer than min_genes core genes"}
    return mito, why


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--assembly", required=True); ap.add_argument("--tblastn", default=None)
    ap.add_argument("--out-mito", required=True); ap.add_argument("--out-nuclear", required=True); ap.add_argument("--json", required=True)
    ap.add_argument("--max-len", type=int, default=250000); ap.add_argument("--min-genes", type=int, default=2); ap.add_argument("--max-gc", type=float, default=0.40)
    a = ap.parse_args()
    seqs = read_fasta(a.assembly)
    mito, why = classify(seqs, genes_per_contig(a.tblastn), a.max_len, a.min_genes, a.max_gc)
    with open(a.out_mito, "w") as fm, open(a.out_nuclear, "w") as fn:
        for name, seq in seqs.items():
            (fm if name in mito else fn).write(f">{name}\n")
            for i in range(0, len(seq), 80):
                (fm if name in mito else fn).write(seq[i:i + 80] + "\n")
    doc = {"sample": a.sample, "n_contigs": len(seqs), "mito_contigs": mito, "mito_bp": sum(len(seqs[c]) for c in mito),
           "nuclear_contigs": len(seqs) - len(mito), "candidates": why,
           "rule": f"core genes >= {a.min_genes}, length <= {a.max_len}, GC < {a.max_gc}"}
    json.dump(doc, open(a.json, "w"), indent=2)
    print(f"[mito_extract] {a.sample}: {len(mito)} mitochondrial contig(s), {doc['mito_bp']} bp; {doc['nuclear_contigs']} nuclear")


if __name__ == "__main__":
    main()
