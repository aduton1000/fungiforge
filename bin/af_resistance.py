#!/usr/bin/env python3
"""FungiForge antifungal-resistance caller (bespoke — there is no ResFinder for fungi).

For the species-relevant target genes in the curated panel it:
  * locates the isolate ortholog among the Funannotate proteins (k-mer prefilter -> pairwise align),
  * reads the residues at the panel's hotspot / known-mutation positions in the reference frame,
  * classifies each as a KNOWN resistance mutation, a novel hotspot variant, or wild-type,
  * handles loss-of-function (FUR1/FCY1/ERG3…), gain-of-function regulators (TAC1/MRR1/UPC2/PDR1),
    efflux over-expression (report presence; needs expression data), and
  * runs the dedicated A. fumigatus cyp51A promoter TR34/TR46 detector (structural).

Every call carries a confidence flag: ONT-only-polished assemblies are marked *provisional*
because homopolymer indels can mimic frameshifts/substitutions exactly where these mutations
live; hybrid- and Illumina-only assemblies are *high* confidence (short-read base accuracy has
no homopolymer-indel problem). Degrades gracefully: without a reference sequence (FungAMR not yet staged) it reports the
gene as searched-but-unresolved rather than guessing, and the structural TR scan still runs.
"""
from __future__ import annotations
import argparse
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cyp51a_TR


# ---------- IO helpers ----------
def read_fasta(path):
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


def read_panel(path):
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts; continue
            rows.append(dict(zip(header, parts)))
    return rows


def read_species(path):
    try:
        first = open(path).readline().strip()
    except Exception:
        return "unknown"
    return first.split("\t")[0].strip() or "unknown"


# ---------- reference index (FungAMR / bundled) ----------
def load_reference_index(data_dir, ref_faa, bundled):
    """gene(lower) -> list of (organism_hint, header, seq)."""
    idx = {}
    faas = []
    if ref_faa and os.path.exists(ref_faa):
        faas.append(ref_faa)
    if data_dir:
        faas += glob.glob(os.path.join(data_dir, "fungamr", "**", "*.f*a*"), recursive=True)
    if bundled and os.path.exists(bundled):
        faas.append(bundled)
    gene_kw = re.compile(r"(cyp51[ab]?|erg1?1|erg[236]|fks[12]|fur1|fcy[12]|hmg1|tac1|mrr1|upc2|pdr1|cdr[12]|mdr1|atrf)", re.I)
    for fa in faas:
        try:
            for h, s in read_fasta(fa).items():
                # a FungAMR reference header can carry several gene names (e.g.
                # "Cdr1_Erg11_Fcy1_Fks1__ACC__Species") — index under EVERY gene keyword,
                # not just the first, or the ERG11 reference would be missed on an ERG11 query.
                genes = set()
                for m in gene_kw.finditer(h):
                    g = m.group(1).lower()
                    if g == "cyp51":
                        g = "cyp51a" if "cyp51a" in h.lower() else "cyp51"
                    genes.add(g)
                for g in genes:
                    idx.setdefault(g, []).append((h, h, s))
        except Exception:
            continue
    return idx


# ---------- alignment / position mapping ----------
def _kmers(s, k=8):
    return {s[i:i + k] for i in range(0, max(0, len(s) - k + 1))}


def best_ortholog(refseq, proteins, top=3):
    """k-mer prefilter then pairwise-align top candidates; return (name, seq, pid) or None."""
    rk = _kmers(refseq)
    if not rk:
        return None
    scored = sorted(((len(rk & _kmers(s)), n) for n, s in proteins.items()), reverse=True)[:top]
    try:
        from Bio.Align import PairwiseAligner
    except Exception:
        return None
    aligner = PairwiseAligner(); aligner.mode = "global"
    aligner.open_gap_score = -10; aligner.extend_gap_score = -0.5
    try:
        from Bio.Align import substitution_matrices
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception:
        aligner.match_score = 2; aligner.mismatch_score = -1
    best = None
    for _, n in scored:
        q = proteins[n]
        try:
            aln = aligner.align(refseq, q)[0]
        except Exception:
            continue
        ident = _identity(aln)
        if best is None or ident > best[2]:
            best = (n, q, ident, aln)
    if best and best[2] >= 0.55:
        return best
    return None


def _identity(aln):
    t, q = aln.aligned
    if len(t) == 0:
        return 0.0
    match = tot = 0
    A, B = aln.target, aln.query
    for (ts, te), (qs, qe) in zip(t, q):
        for i in range(te - ts):
            tot += 1
            if A[ts + i] == B[qs + i]:
                match += 1
    return match / tot if tot else 0.0


def ref_to_query_residue(aln, ref_pos1):
    """Residue in the query aligned to 1-based reference position ref_pos1, or None (deletion)."""
    p = ref_pos1 - 1
    t, q = aln.aligned
    for (ts, te), (qs, qe) in zip(t, q):
        if ts <= p < te:
            return aln.query[qs + (p - ts)]
    return None


# ---------- calling ----------
def call_substitutions(gene, refseq, aln, row):
    """Yield resistance/variant calls at hotspot + known positions."""
    known = {}
    for tok in (row.get("known_mutations", "") or "").split(","):
        tok = tok.strip()
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", tok)
        if m:
            known.setdefault(int(m.group(2)), []).append((m.group(1), m.group(3)))
    positions = set(known)
    for h in (row.get("hotspot_aa", "") or "").split(","):
        h = h.strip()
        if h.isdigit():
            positions.add(int(h))
    calls = []
    for pos in sorted(positions):
        if pos > len(refseq):
            continue
        wt = refseq[pos - 1]
        obs = ref_to_query_residue(aln, pos)
        if obs is None or obs == wt or obs == "-":
            continue
        change = f"{wt}{pos}{obs}"
        is_known = any(mut == obs for (_, mut) in known.get(pos, []))
        calls.append({"gene": gene, "change": change, "known": is_known,
                      "class": "known_resistance_mutation" if is_known else "novel_hotspot_variant"})
    return calls


def main():
    ap = argparse.ArgumentParser(description="FungiForge antifungal-resistance caller")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--assembly")
    ap.add_argument("--gbk")
    ap.add_argument("--panel", required=True)
    ap.add_argument("--data-dir", default="")
    ap.add_argument("--ref-faa", help="explicit reference protein FASTA (headers containing gene names)")
    ap.add_argument("--polish-mode", default="unknown",
                    choices=["hybrid", "illumina_only", "ont_only", "unknown"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    species = read_species(a.species)
    panel = read_panel(a.panel)
    proteins = read_fasta(a.proteins) if os.path.exists(a.proteins) else {}
    bundled_ref = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                               "fungiforge", "resources", "af_reference_proteins.faa")
    refidx = load_reference_index(a.data_dir, a.ref_faa, bundled_ref)

    # Illumina base accuracy carries no homopolymer-indel risk, so short-read-only
    # assemblies are high-confidence for point-mutation resistance calls, like hybrid.
    conf = {"hybrid": "high", "illumina_only": "high",
            "ont_only": "provisional_ont_only", "unknown": "provisional"}[a.polish_mode]
    # species relevance: organism_regex matches the species, or genus matches, or 'spp.'
    genus = species.split()[0] if species and species != "unknown" else ""
    def relevant(row):
        org = row.get("organism_regex", "")
        return bool(genus) and (re.search(re.escape(species), org, re.I) or
                                re.search(re.escape(genus), org, re.I) or "spp." in org)

    calls, searched = [], []
    for row in panel:
        gene = row["gene"]
        mech = row.get("mechanism", "")
        if not relevant(row):
            continue
        searched.append(gene)
        if mech == "promoter_TR":
            continue  # handled below via cyp51a_TR
        refs = refidx.get(gene.lower()) or refidx.get(re.sub(r"[_].*$", "", gene.lower()))
        if not refs:
            calls.append({"gene": gene, "drug_class": row.get("drug_class"), "mechanism": mech,
                          "status": "no_reference", "confidence": conf,
                          "note": "gene searched; no reference sequence staged (install FungAMR) — cannot resolve residues"})
            continue
        # Pick a reference from the isolate's SPECIES (headers are GENE__ACC__Species). The
        # panel's hotspot numbering is species-specific, so a wrong-species reference would
        # misnumber every residue. Fall back to genus, then to the longest available.
        toks = [t for t in re.split(r"\s+", species.lower()) if len(t) > 2]
        sp_refs = [s for (h, _, s) in refs if toks and all(t in h.lower() for t in toks)]
        if not sp_refs and toks:
            sp_refs = [s for (h, _, s) in refs if toks[0] in h.lower()]
        refseq = max(sp_refs, key=len) if sp_refs else max((s for _, _, s in refs), key=len)
        orth = best_ortholog(refseq, proteins) if proteins else None
        if not orth:
            calls.append({"gene": gene, "drug_class": row.get("drug_class"), "mechanism": mech,
                          "status": "not_detected", "confidence": conf})
            continue
        name, qseq, pid, aln = orth
        base = {"gene": gene, "drug_class": row.get("drug_class"), "drugs": row.get("drugs"),
                "mechanism": mech, "ortholog": name, "identity": round(pid, 3),
                "confidence": conf, "source": row.get("source", "")}
        if mech == "substitution":
            subs = call_substitutions(gene, refseq, aln, row)
            if subs:
                for c in subs:
                    calls.append({**base, "status": "variant", **c})
            else:
                calls.append({**base, "status": "wild_type"})
        elif mech == "loss_of_function":
            trunc = len(qseq) < 0.85 * len(refseq)
            calls.append({**base, "status": "loss_of_function" if trunc else "intact",
                          "note": f"length {len(qseq)}/{len(refseq)} aa" + (" — truncated (candidate LoF)" if trunc else "")})
        elif mech == "GOF":
            subs = call_substitutions(gene, refseq, aln, row)
            calls.append({**base, "status": "present",
                          "gof_variants": [c["change"] for c in subs if c["known"]],
                          "note": "regulator; gain-of-function inferred from known residues (confirm with expression/phenotype)"})
        elif mech == "overexpression":
            calls.append({**base, "status": "present",
                          "note": "efflux/target present; over-expression needs RNA-seq — not callable from genome alone"})

    # cyp51A promoter tandem repeats (A. fumigatus)
    tr = None
    if re.search(r"Aspergillus\s+fumigatus", species, re.I) and a.assembly:
        tr = cyp51a_TR.detect_tr(a.assembly, a.gbk)
        if tr.get("tr_detected"):
            calls.append({"gene": "cyp51A_promoter", "drug_class": "azole", "mechanism": "promoter_TR",
                          "status": "resistance", "change": tr["tr_type"], "known": True,
                          "confidence": "high" if a.polish_mode in ("hybrid", "illumina_only") else "medium",
                          "note": tr["note"] + " — pairs with cyp51A L98H (TR34) or Y121F+T289A (TR46)"})

    resistant_classes = sorted({c.get("drug_class") for c in calls
                                if c.get("known") or c.get("status") in ("loss_of_function", "resistance")
                                and c.get("drug_class")})
    out = {"sample": a.sample, "stage": "resistance", "species": species,
           "polish_mode": a.polish_mode, "genes_searched": sorted(set(searched)),
           "calls": calls, "cyp51A_TR": tr,
           "summary": {"resistant_drug_classes": [c for c in resistant_classes if c],
                       "n_known_mutations": sum(1 for c in calls if c.get("known")),
                       "reference_available": bool(refidx)}}
    json.dump(out, open(a.out, "w"), indent=2)
    print(f"[af_resistance] {a.sample} · {species} · {len(calls)} calls · "
          f"known={out['summary']['n_known_mutations']} · ref={'yes' if refidx else 'no'}")


if __name__ == "__main__":
    main()
