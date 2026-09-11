#!/usr/bin/env python3
"""Detect Aspergillus fumigatus cyp51A promoter tandem repeats (TR34 / TR46 / TR46^3 / TR53).

These are the dominant *environmental* azole-resistance mechanism and are structural
(promoter insertions), not point mutations — so they are invisible to a substitution
caller and need a dedicated structural scan. Approach:

  1. locate cyp51A in the annotated GBK (gene/product name), get its contig + start + strand;
  2. extract the ~600 bp upstream promoter from the assembly (strand-aware);
  3. scan the promoter for a tandem duplication of a ~34/46/53 bp unit (the region between
     the two SrbA-binding sites). Report the unit length, copy number and TR type.

The canonical genotype pairings (annotated by af_resistance.py) are TR34/L98H and
TR46/Y121F+T289A. Usable as a module (`detect_tr`) or standalone.
"""
from __future__ import annotations
import argparse
import json
import re

CYP51A_NAMES = re.compile(r"cyp51a|cyp51\b|erg11|eburicol|sterol.{0,3}14.?alpha.?demethyl", re.I)
KNOWN_TR = {34: "TR34", 46: "TR46", 53: "TR53"}


def _read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name: seqs[name] = "".join(buf)
    return seqs


def _revcomp(s):
    return s.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


CYP51A_GENE = re.compile(r"^(cyp51a|erg11|afua_4g06890|afu4g06890)$", re.I)


def _locate_cyp51a(gbk_path):
    """Return (contig, start, end, strand, contig_seq) of cyp51A from a GenBank file, or None.

    Two passes: an explicit gene name (CYP51A / erg11 / the A. fumigatus locus tag) wins;
    only if none is present fall back to product/note wording, which also matches the
    CYP51B paralog ("Sterol 14-alpha demethylase") and must not be preferred.
    The record's own sequence is returned because funannotate renames/sorts contigs, so
    GBK contig names do not map onto the assembly FASTA the pipeline holds.
    """
    try:
        from Bio import SeqIO
    except Exception:
        return None
    fallback = None
    for rec in SeqIO.parse(gbk_path, "genbank"):
        for feat in rec.features:
            if feat.type not in ("gene", "mRNA", "CDS"):
                continue
            q = feat.qualifiers
            hit = (rec.id, int(feat.location.start), int(feat.location.end),
                   1 if feat.location.strand in (1, None) else -1, str(rec.seq))
            if any(CYP51A_GENE.match(g.strip()) for g in q.get("gene", []) + q.get("locus_tag", [])):
                return hit
            names = " ".join(q.get("gene", []) + q.get("product", []) + q.get("note", []))
            if fallback is None and CYP51A_NAMES.search(names) and not re.search(r"cyp51b", names, re.I):
                fallback = hit
    return fallback


def _tandem_scan(promoter, unit_lengths=range(20, 61), min_id=0.90):
    """Find the best tandem array in the promoter. Returns (unit_len, copies, pos) or None."""
    best = None
    n = len(promoter)
    for L in unit_lengths:
        i = 0
        while i + 2 * L <= n:
            unit = promoter[i:i + L]
            copies, j = 1, i + L
            while j + L <= n:
                nxt = promoter[j:j + L]
                ident = sum(a == b for a, b in zip(unit, nxt)) / L
                if ident >= min_id:
                    copies += 1; j += L
                else:
                    break
            if copies >= 2:
                span = copies * L
                if best is None or span > best[1] * best[0]:
                    best = (L, copies, i)
                i = j
            else:
                i += 1
    return best


def detect_tr(assembly_path, gbk_path=None, promoter_bp=600):
    out = {"tr_detected": False, "tr_type": "none", "unit_len": None,
           "copies": None, "cyp51A_located": False, "note": ""}
    loc = _locate_cyp51a(gbk_path) if gbk_path else None
    contigs = _read_fasta(assembly_path)
    if not loc:
        out["note"] = "cyp51A not located in GBK; TR scan skipped (provide annotated GBK)"
        return out
    out["cyp51A_located"] = True
    contig, start, end, strand, gbk_seq = loc
    # Coordinates are in GBK space -> use the GBK record's sequence; the assembly FASTA is
    # only a fallback (its contig names are Flye's, funannotate's are size-sorted).
    seq = gbk_seq or contigs.get(contig, "")
    if not seq:
        out["note"] = f"contig {contig} has no sequence in GBK and was not found in assembly"
        return out
    if strand == 1:
        promoter = seq[max(0, start - promoter_bp):start]
    else:
        promoter = _revcomp(seq[end:end + promoter_bp])
    hit = _tandem_scan(promoter.upper())
    if hit:
        L, copies, pos = hit
        # snap to the nearest canonical TR unit length for naming
        canon = min(KNOWN_TR, key=lambda k: abs(k - L))
        tr_type = KNOWN_TR[canon]
        if copies >= 3 and canon == 46:
            tr_type = "TR46_3"
        out.update(tr_detected=True, tr_type=tr_type, unit_len=L, copies=copies,
                   note=f"tandem array of {copies}×{L} bp in cyp51A promoter (~{tr_type})")
    else:
        out["note"] = "no tandem repeat detected in cyp51A promoter (wild-type promoter)"
    return out


def main():
    ap = argparse.ArgumentParser(description="A. fumigatus cyp51A promoter TR detector")
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--gbk")
    ap.add_argument("--out")
    a = ap.parse_args()
    res = detect_tr(a.assembly, a.gbk)
    if a.out:
        json.dump(res, open(a.out, "w"), indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
