#!/usr/bin/env python3
"""Extract the secondary identification loci from an assembly (W2.3, stage 08).

Protein-coding markers (CaM, BenA, TEF1, RPB2) are located with tblastn hits of the bundled
reference proteins (fungiforge/resources/markers/marker_reference_proteins.faa); the LSU D1/D2
region comes from barrnap's 28S coordinates. Each locus is written in coding orientation with
flanks, ready for a nucleotide search against the type-material reference set of that locus.

  extract_markers.py --sample ID --assembly asm.fa --tblastn tblastn.tsv [--rrna-gff rrna.gff]
                     [--flank 300] [--max-gap 5000] [--min-bitscore 80] [--lsu-len 900]
                     --out-prefix ID.marker --json ID.markers.json

tblastn.tsv: outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
qseqid = "<MARKER>|<accession>|<organism>". HSPs of the best-scoring contig/strand are chained
when they lie within --max-gap of each other (introns), and the chained span ± flank is the locus.
"""
from __future__ import annotations
import argparse, json, os

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


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


def parse_tblastn(path):
    """-> {marker: [hsp dicts]} with 1-based subject coordinates and strand."""
    hits = {}
    if not path or not os.path.exists(path):
        return hits
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            marker = f[0].split("|")[0]
            ss, se = int(f[6]), int(f[7])
            strand = "+" if se >= ss else "-"
            hits.setdefault(marker, []).append({
                "contig": f[1], "pident": float(f[2]), "alen": int(f[3]), "qstart": int(f[4]), "qend": int(f[5]),
                "sstart": min(ss, se), "send": max(ss, se), "strand": strand, "evalue": float(f[8]), "bitscore": float(f[9])})
    return hits


def chain_best(hsps, max_gap=5000, min_bitscore=80.0):
    """Best contig/strand by summed bitscore; chain HSPs within max_gap of the top HSP; return locus."""
    hsps = [h for h in hsps if h["bitscore"] >= min_bitscore]
    if not hsps:
        return None
    groups = {}
    for h in hsps:
        groups.setdefault((h["contig"], h["strand"]), []).append(h)
    key = max(groups, key=lambda k: sum(h["bitscore"] for h in groups[k]))
    g = sorted(groups[key], key=lambda h: -h["bitscore"])
    chain = [g[0]]
    changed = True
    while changed:
        changed = False
        lo, hi = min(h["sstart"] for h in chain), max(h["send"] for h in chain)
        for h in g:
            if h in chain:
                continue
            if h["send"] >= lo - max_gap and h["sstart"] <= hi + max_gap:
                chain.append(h); changed = True
    qcov = len({q for h in chain for q in range(h["qstart"], h["qend"] + 1)})
    return {"contig": key[0], "strand": key[1], "start": min(h["sstart"] for h in chain), "end": max(h["send"] for h in chain),
            "bitscore": round(sum(h["bitscore"] for h in chain), 1), "n_hsps": len(chain),
            "pident_best": max(h["pident"] for h in chain), "query_residues_covered": qcov,
            "n_contigs_hit": len({c for c, _ in groups})}


def lsu_from_gff(gff, contigs, lsu_len=900):
    """Longest 28S feature from a barrnap GFF -> D1/D2 region (first lsu_len bp of the 5' end)."""
    best = None
    if not gff or not os.path.exists(gff):
        return None
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or "28S" not in f[8]:
                continue
            c, s, e, strand = f[0], int(f[3]), int(f[4]), f[6]
            if c not in contigs:
                continue
            if best is None or (e - s) > (best["end"] - best["start"]):
                best = {"contig": c, "start": s, "end": e, "strand": strand, "partial": "partial" in f[8].lower()}
    if not best:
        return None
    seq = contigs[best["contig"]][best["start"] - 1:best["end"]]
    if best["strand"] == "-":
        seq = revcomp(seq)
    best["full_len"] = len(seq)
    best["seq"] = seq[:lsu_len]
    best["truncated"] = len(seq) < lsu_len
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--assembly", required=True)
    ap.add_argument("--tblastn", default=None); ap.add_argument("--rrna-gff", default=None)
    ap.add_argument("--flank", type=int, default=300); ap.add_argument("--max-gap", type=int, default=5000)
    ap.add_argument("--min-bitscore", type=float, default=80.0); ap.add_argument("--lsu-len", type=int, default=900)
    ap.add_argument("--markers", default="CaM,BenA,TEF1,RPB2")
    ap.add_argument("--out-prefix", required=True); ap.add_argument("--json", required=True)
    a = ap.parse_args()
    contigs = read_fasta(a.assembly)
    hits = parse_tblastn(a.tblastn)
    summary = {"sample": a.sample, "loci": {}}
    for m in a.markers.split(","):
        loc = chain_best(hits.get(m, []), a.max_gap, a.min_bitscore)
        out = f"{a.out_prefix}.{m}.fasta"
        if not loc or loc["contig"] not in contigs:
            summary["loci"][m] = {"found": False}
            open(out, "w").close()
            continue
        s = max(1, loc["start"] - a.flank); e = min(len(contigs[loc["contig"]]), loc["end"] + a.flank)
        seq = contigs[loc["contig"]][s - 1:e]
        if loc["strand"] == "-":
            seq = revcomp(seq)
        with open(out, "w") as fh:
            fh.write(f">{a.sample}|{m}|{loc['contig']}:{s}-{e}({loc['strand']}) bitscore={loc['bitscore']} hsps={loc['n_hsps']}\n{seq}\n")
        summary["loci"][m] = {"found": True, "file": os.path.basename(out), "extracted_start": s, "extracted_end": e, "length": len(seq), **loc}
    lsu = lsu_from_gff(a.rrna_gff, contigs, a.lsu_len)
    out = f"{a.out_prefix}.LSU.fasta"
    if lsu:
        with open(out, "w") as fh:
            fh.write(f">{a.sample}|LSU|{lsu['contig']}:{lsu['start']}-{lsu['end']}({lsu['strand']}) d1d2_len={len(lsu['seq'])} full_28S_len={lsu['full_len']}\n{lsu['seq']}\n")
        summary["loci"]["LSU"] = {"found": True, "file": os.path.basename(out), "contig": lsu["contig"], "start": lsu["start"],
                                  "end": lsu["end"], "strand": lsu["strand"], "length": len(lsu["seq"]), "full_28S_len": lsu["full_len"],
                                  "truncated": lsu["truncated"]}
    else:
        open(out, "w").close()
        summary["loci"]["LSU"] = {"found": False}
    with open(a.json, "w") as fh:
        json.dump(summary, fh, indent=2)
    found = [m for m, v in summary["loci"].items() if v["found"]]
    print(f"[extract_markers] {a.sample}: found {', '.join(found) or 'none'}; missing {', '.join(m for m in summary['loci'] if m not in found) or 'none'}")


if __name__ == "__main__":
    main()
