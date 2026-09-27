#!/usr/bin/env python3
"""Annotate and characterise the mitochondrial genome (W2.7, stage 04b).

  mito_annotate.py --sample ID --mito ID.mito.fasta --tblastn genes.tsv [--rrna-blastn rrna.tsv]
                   [--mpileup mito.pileup] [--mito-depth D] [--nuclear-depth D]
                   --out-json ID.organelle.json --out-gff ID.mito.gff

genes.tsv: tblastn of the bundled core proteins vs the mito contigs (outfmt 6 + qlen); rrna.tsv:
blastn of the Af293 rRNAs (outfmt 6). Per gene the best species' HSPs on the best contig/strand
are chained (gaps >= --min-intron bp between consecutive HSPs = introns), giving presence, copies
(distinct loci), intron count and coverage of the reference protein. Circularity = the contig's
end repeats its start (>= --overlap bp). Heteroplasmy = pileup positions (depth >= 20) where a
minor allele reaches >= 10 %. Copy number = mito depth / nuclear depth."""
from __future__ import annotations
import argparse, collections, json, os

CORE = ["cox1", "cox2", "cox3", "cob", "atp6", "atp8", "atp9", "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6", "rps3"]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    if not path or not os.path.exists(path):
        return seqs
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


def parse_hits(path, with_qlen=True):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            try:
                ss, se = int(f[6]), int(f[7])
                rows.append({"query": f[0], "gene": f[0].split("|")[0], "species": f[0].split("|")[-1], "contig": f[1], "pident": float(f[2]),
                             "alen": int(f[3]), "qstart": int(f[4]), "qend": int(f[5]), "sstart": min(ss, se), "send": max(ss, se),
                             "strand": "+" if se >= ss else "-", "evalue": float(f[8]), "bitscore": float(f[9]),
                             "qlen": int(f[10]) if with_qlen and len(f) > 10 else None})
            except ValueError:
                continue
    return rows


def chain(hsps, max_gap=6000):
    """Group HSPs of one query on one contig/strand into loci; return list of loci dicts."""
    hsps = sorted(hsps, key=lambda h: h["sstart"])
    loci = []
    for h in hsps:
        if loci and h["sstart"] - loci[-1]["end"] <= max_gap:
            L = loci[-1]
            L["end"] = max(L["end"], h["send"]); L["hsps"].append(h); L["bitscore"] += h["bitscore"]
        else:
            loci.append({"start": h["sstart"], "end": h["send"], "hsps": [h], "bitscore": h["bitscore"]})
    return loci


def annotate_genes(rows, min_intron=100, max_evalue=1e-5):
    """Per core gene: best locus (species with the highest summed bitscore), copies, introns, coverage."""
    out = {}
    rows = [r for r in rows if r["evalue"] <= max_evalue]
    for gene in CORE:
        g = [r for r in rows if r["gene"] == gene]
        if not g:
            out[gene] = {"present": False}
            continue
        best = None
        by_query = collections.defaultdict(list)
        for r in g:
            by_query[(r["query"], r["contig"], r["strand"])].append(r)
        for (q, contig, strand), hs in by_query.items():
            for L in chain(hs):
                if best is None or L["bitscore"] > best["bitscore"]:
                    best = dict(L, query=q, contig=contig, strand=strand)
        hs = sorted(best["hsps"], key=lambda h: h["sstart"])
        introns = [hs[i + 1]["sstart"] - hs[i]["send"] - 1 for i in range(len(hs) - 1)]
        introns = [x for x in introns if x >= min_intron]
        qlen = hs[0]["qlen"]
        covered = len({p for h in hs for p in range(h["qstart"], h["qend"] + 1)})
        # copies = loci of this gene on ONE contig (a real duplication); contigs_with_gene counts the
        # contigs carrying it (> 1 when the assembler left overlapping fragments of the circle)
        per_contig = collections.Counter()
        for (q, contig, strand), hs2 in by_query.items():
            if q != best["query"]:
                continue
            per_contig[contig] += len(chain(hs2))
        out[gene] = {"present": True, "contig": best["contig"], "start": best["start"], "end": best["end"], "strand": best["strand"],
                     "reference": best["query"].split("|")[-1], "pident": round(max(h["pident"] for h in hs), 1),
                     "query_coverage": round(covered / qlen, 2) if qlen else None, "n_exons": len(hs), "n_introns": len(introns),
                     "intron_lengths": introns, "copies": max(per_contig.values()) if per_contig else 1, "contigs_with_gene": len(per_contig)}
    return out


def annotate_rrna(rows, min_len=200):
    out = {}
    for name in ("rnl", "rns"):
        g = [r for r in rows if r["gene"] == name and r["alen"] >= min_len]
        if not g:
            out[name] = {"present": False}; continue
        best = max(g, key=lambda r: r["bitscore"])
        out[name] = {"present": True, "contig": best["contig"], "start": best["sstart"], "end": best["send"], "strand": best["strand"], "pident": round(best["pident"], 1), "alen": best["alen"]}
    return out


def circular(seq, overlap=50, max_overlap=2000):
    """True when the contig end repeats its start (an unresolved circular junction)."""
    s = seq.upper()
    for k in range(min(max_overlap, len(s) // 2), overlap - 1, -1):
        if s[:k] == s[-k:]:
            return k
    return 0


def heteroplasmy(pileup_path, min_depth=20, min_frac=0.10):
    """samtools mpileup (no -s): chrom pos ref depth bases quals -> sites where a minor SUBSTITUTION
    allele reaches >= min_frac (indels are ignored: ONT homopolymer errors would otherwise dominate;
    use min_frac 0.20 for ONT, 0.10 for Illumina)."""
    sites, n_cov = [], 0
    if not pileup_path or not os.path.exists(pileup_path):
        return None
    with open(pileup_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                depth = int(f[3])
            except ValueError:
                continue
            if depth < min_depth:
                continue
            n_cov += 1
            bases = f[4]
            counts = collections.Counter()
            i = 0
            while i < len(bases):
                c = bases[i]
                if c == "^":
                    i += 2; continue
                if c == "$":
                    i += 1; continue
                if c in "+-":
                    j = i + 1; num = ""
                    while j < len(bases) and bases[j].isdigit():
                        num += bases[j]; j += 1
                    i = j + int(num or 0); counts["indel"] += 1; continue
                if c in ".,":
                    counts["ref"] += 1
                elif c.upper() in "ACGT":
                    counts[c.upper()] += 1
                i += 1
            counts.pop("indel", None)
            tot = sum(counts.values())
            if tot < min_depth:
                continue
            top2 = counts.most_common(2)
            if len(top2) > 1 and top2[1][1] / tot >= min_frac:
                sites.append({"contig": f[0], "pos": int(f[1]), "depth": tot, "alleles": {k: round(v / tot, 3) for k, v in top2}})
    return {"n_positions_covered": n_cov, "n_heteroplasmic_sites": len(sites), "min_minor_frac": min_frac, "sites": sites[:200]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--mito", required=True); ap.add_argument("--tblastn", default=None)
    ap.add_argument("--rrna-blastn", default=None); ap.add_argument("--mpileup", default=None)
    ap.add_argument("--mito-depth", type=float, default=None); ap.add_argument("--nuclear-depth", type=float, default=None)
    ap.add_argument("--min-intron", type=int, default=100); ap.add_argument("--het-min-frac", type=float, default=0.10)
    ap.add_argument("--out-json", required=True); ap.add_argument("--out-gff", required=True)
    a = ap.parse_args()
    seqs = read_fasta(a.mito)
    size = sum(len(s) for s in seqs.values())
    doc = {"sample": a.sample, "stage": "organelle", "mito_present": bool(seqs), "n_contigs": len(seqs), "mito_size": size,
           "contigs": {n: {"length": len(s), "gc": round(sum(1 for c in s.upper() if c in "GC") / max(1, len(s)), 3), "circular_overlap": circular(s)} for n, s in seqs.items()}}
    if not seqs:
        doc.update(genes={}, rrna={}, n_core_genes=0, note="no mitochondrial contig separated by stage 04")
        json.dump(doc, open(a.out_json, "w"), indent=2); open(a.out_gff, "w").write("##gff-version 3\n")
        print(f"[mito_annotate] {a.sample}: no mitochondrial contigs"); return
    genes = annotate_genes(parse_hits(a.tblastn), a.min_intron)
    rrna = annotate_rrna(parse_hits(a.rrna_blastn, with_qlen=False))
    doc["genes"] = genes; doc["rrna"] = rrna
    doc["n_core_genes"] = sum(1 for g in genes.values() if g.get("present"))
    doc["missing_core_genes"] = [g for g in CORE if not genes[g].get("present")]
    doc["n_introns"] = sum(g.get("n_introns", 0) for g in genes.values() if g.get("present"))
    doc["genes_with_extra_copies"] = [g for g, v in genes.items() if v.get("present") and v.get("copies", 1) > 1]
    # overlapping fragments of the circle (Flye often leaves 2-3): most genes sit on more than one contig
    present = [v for v in genes.values() if v.get("present")]
    multi = sum(1 for v in present if v.get("contigs_with_gene", 1) > 1)
    doc["redundant_contigs"] = bool(len(seqs) > 1 and present and multi >= 0.5 * len(present))
    longest = max(seqs, key=lambda n: len(seqs[n]))
    doc["longest_contig"] = {"name": longest, "length": len(seqs[longest])}
    doc["mito_size_estimate"] = len(seqs[longest]) if doc["redundant_contigs"] else size
    if doc["redundant_contigs"]:
        doc["note"] = "the mitochondrial contigs overlap (most core genes occur on more than one contig): the assembler left the circular genome as redundant fragments; size estimate = longest contig"
    doc["circular"] = any(c["circular_overlap"] >= 50 for c in doc["contigs"].values()) and len(seqs) == 1
    doc["heteroplasmy"] = heteroplasmy(a.mpileup, min_frac=a.het_min_frac)
    doc["mito_depth"] = a.mito_depth; doc["nuclear_depth"] = a.nuclear_depth
    doc["copy_ratio"] = round(a.mito_depth / a.nuclear_depth, 1) if a.mito_depth and a.nuclear_depth else None
    with open(a.out_gff, "w") as fh:
        fh.write("##gff-version 3\n")
        for name, g in genes.items():
            if g.get("present"):
                fh.write(f"{g['contig']}\tfungiforge\tgene\t{g['start']}\t{g['end']}\t{g['pident']}\t{g['strand']}\t.\tID={name};Name={name};introns={g['n_introns']};ref={g['reference']}\n")
        for name, r in rrna.items():
            if r.get("present"):
                fh.write(f"{r['contig']}\tfungiforge\trRNA\t{r['start']}\t{r['end']}\t{r['pident']}\t{r['strand']}\t.\tID={name};Name={name}\n")
    json.dump(doc, open(a.out_json, "w"), indent=2)
    print(f"[mito_annotate] {a.sample}: {size} bp in {len(seqs)} contig(s) (estimate {doc['mito_size_estimate']}{', redundant' if doc['redundant_contigs'] else ''}), {doc['n_core_genes']}/15 core genes, {doc['n_introns']} introns, "
          f"circular={doc['circular']}, heteroplasmic sites={(doc['heteroplasmy'] or {}).get('n_heteroplasmic_sites', 'NA')}, copy ratio={doc['copy_ratio']}")


if __name__ == "__main__":
    main()
