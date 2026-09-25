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

W2.4 additions: (1) the panel is the curated overlay PLUS rows derived from the FungAMR
catalogue (bin/fungamr_panel.py) so every known substitution carries an evidence tier
(1 strongest .. 8 = seen in a resistant natural isolate only; tier-8-only changes are reported
as `associated_unvalidated`, not as resistance); (2) with `--bam` (reads mapped to the GenBank
records, bin/gbk_to_fasta.py) every hotspot is re-genotyped from the reads
(bin/read_genotype.py): allele frequency, zygosity, agreement with the assembly, plus known
alleles the assembly missed (heterozygous in a diploid), the cyp51A TR site as a read-level
insertion, and locus/genome depth ratio as a copy-number signal for target and efflux genes;
(3) confidence combines polish mode, read support and evidence tier.
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
import fungamr_panel
import read_genotype


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
def panel_positions(row):
    """{position: [(wt, mut), ...]} of known mutations and the set of hotspot positions of a row."""
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
    return known, positions


def evidence_for(evidence, gene, species, change):
    """FungAMR tier/classes for a change, trying the species, then any species of the genus.
    The returned dict carries `species`: the organism the entry was recorded in. Residue numbering
    is species-specific, so an entry from ANOTHER species of the genus is association at best and
    call_substitutions never counts it as a known mutation."""
    key = (gene.lower(), (species or "").lower(), change)
    if key in evidence:
        return dict(evidence[key], species=species)
    genus = (species or "").split()[0].lower() if species else ""
    for (g, sp, ch), v in evidence.items():
        if g == gene.lower() and ch == change and sp.split()[0] == genus:
            return dict(v, species=sp)
    return None


def call_substitutions(gene, refseq, aln, row, evidence=None, species=None, curated=True):
    """Resistance/variant calls at hotspot + known positions of one panel row."""
    known, positions = panel_positions(row)
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
        ev = evidence_for(evidence or {}, gene, species, change) if is_known else None
        curated_known = row.get("curated_known")            # set by merge_panels: the overlay's own mutations
        in_overlay = curated and (change in curated_known if curated_known is not None else True)
        tier = "curated" if (is_known and in_overlay) else (ev["tier"] if ev else None)
        c = {"gene": gene, "change": change, "known": is_known,
             "class": "known_resistance_mutation" if is_known else "novel_hotspot_variant",
             "evidence_tier": tier, "ref_pos": pos, "wt": wt, "obs": obs}
        if is_known and tier == 8:
            c["class"], c["known"] = "associated_unvalidated", False
            c["note"] = "FungAMR tier 8: seen in a resistant isolate without validation — association only"
        if ev and species and ev.get("species") and ev["species"].lower() != species.lower() and tier != "curated":
            # the only evidence is an entry recorded in another species of the genus: the residue
            # numbering is that species', so this is at most an association, never a known mutation
            c["class"], c["known"] = "associated_other_species", False
            c["evidence_species"] = ev["species"]
            c["note"] = f"evidence recorded in {ev['species']}, not in this species — association only"
        if ev and ev.get("classes"):
            c["evidence_classes"] = ev["classes"]
        calls.append(c)
    return calls


def merge_panels(curated_rows, fungamr_rows):
    """Curated overlay first; FungAMR rows add positions/mutations to a matching (gene, organism)
    row or become new substitution rows. Returns rows with a `curated` flag."""
    out = []
    index = {}
    for r in curated_rows:
        r = dict(r, curated=True, curated_known={t for t in (r.get("known_mutations") or "").split(",") if t and t != "NA"})
        out.append(r)
        index[(r["gene"].lower(), (r.get("organism_regex") or "").lower())] = r
    for f in fungamr_rows:
        names = [re.sub(r"\\", "", n).lower() for n in (f.get("organism_regex") or "").split("|") if n]
        target = None
        for r in out:
            if r["gene"].lower() != f["gene"].lower() or r.get("mechanism") != "substitution":
                continue
            org = (r.get("organism_regex") or "").lower()
            if any(n == org or re.search(r"\b" + re.escape(n) + r"\b", org) for n in names):
                target = r; break
        if target:
            km = [t for t in (target.get("known_mutations") or "").split(",") if t and t != "NA"]
            hs = [t for t in (target.get("hotspot_aa") or "").split(",") if t and t != "NA"]
            for t in (f.get("known_mutations") or "").split(","):
                if t and t not in km:
                    km.append(t)
            for t in (f.get("hotspot_aa") or "").split(","):
                if t and t not in hs:
                    hs.append(t)
            target["known_mutations"], target["hotspot_aa"] = ",".join(km), ",".join(hs)
            target["fungamr_merged"] = True
        else:
            out.append(dict(f, curated=False, source=f.get("source", "FungAMR")))
    return out


def combined_confidence(polish_mode, read_support=None, tier=None):
    """high: read-confirmed, or a high-accuracy assembly (hybrid / Illumina); provisional_ont_only:
    ONT-only assembly without read confirmation; discordant: the reads contradict the assembly;
    low_evidence: tier-8 association only (whatever the reads say)."""
    base = {"hybrid": "high", "illumina_only": "high", "ont_only": "provisional_ont_only", "unknown": "provisional"}[polish_mode]
    if read_support:
        if read_support.get("agrees") is False:
            return "discordant"
        if read_support.get("agrees") is True:
            base = "high"
    if tier == 8:
        return "low_evidence"
    return base


def read_support_for(reads, cds, aln, row, refseq, min_reads):
    """Read-level genotype of every panel position of a row, keyed by REFERENCE position.
    The query residue index comes from the protein alignment (reference numbering != isolate
    numbering when the ortholog has indels)."""
    _, positions = panel_positions(row)
    ref_to_query = {}
    t, q = aln.aligned
    for pos in positions:
        p = pos - 1
        for (ts, te), (qs, qe) in zip(t, q):
            if ts <= p < te:
                ref_to_query[pos] = int(qs + (p - ts) + 1)      # plain int: the alignment arrays are numpy
                break
    if not ref_to_query:
        return {}
    try:
        geno = read_genotype.genotype_residues(reads, cds, sorted(set(ref_to_query.values())), min_reads)
    except Exception as e:  # noqa: BLE001
        sys.stderr.write(f"[af_resistance] read genotyping failed for {cds.get('contig')}: {e}\n")
        return {}
    return {pos: dict(geno[qpos], query_pos=qpos) for pos, qpos in ref_to_query.items() if qpos in geno}


def tr_read_check(reads, tr, min_reads, min_frac):
    """Reads spanning the cyp51A TR site, mapped to the assembly: an extra copy in the reads is an
    insertion, a copy the reads lack is a deletion (30-140 bp). call: agrees_with_assembly |
    reads_have_extra_copy | reads_lack_copy | mixed | insufficient."""
    try:
        if isinstance(reads, tuple):
            lines = reads[1]
        else:
            lines = read_genotype.sam_from_bam(reads, tr["contig"], max(1, tr["tr_site_start"] - 300), tr["tr_site_end"] + 300)
        sup = read_genotype.indel_support(lines, tr["contig"], (tr["tr_site_start"], tr["tr_site_end"]), 30, 140)
    except Exception as e:  # noqa: BLE001
        return {"call": "not_run", "error": str(e)[:120]}
    n, ins, dele = sup["spanning"], sup["ins_frac"], sup["del_frac"]
    if n < min_reads or ins is None:
        call = "insufficient"
    elif ins >= min_frac:
        call = "reads_have_extra_copy"
    elif dele >= min_frac:
        call = "reads_lack_copy"
    elif ins <= 1 - min_frac and dele <= 1 - min_frac:
        call = "agrees_with_assembly"
    else:
        call = "mixed"
    return dict(sup, call=call, agrees=(call == "agrees_with_assembly") if call not in ("insufficient", "not_run") else None)


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
    ap.add_argument("--bam", default=None, help="reads mapped to the GenBank records (bin/gbk_to_fasta.py); enables read-level genotyping")
    ap.add_argument("--sam", default=None, help="SAM text instead of --bam (tests)")
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--min-frac", type=float, default=0.8)
    ap.add_argument("--fungamr-min-tier", type=int, default=8, help="keep FungAMR mutations with tier <= this")
    ap.add_argument("--no-fungamr-panel", action="store_true", help="curated overlay only")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    species = read_species(a.species)
    # Panel: the path given (repo copy), else the copy packaged with the `fungiforge`
    # Python package (present inside the container image / cli-env), else fail loudly.
    panel_path = a.panel
    if not os.path.exists(panel_path):
        try:
            from importlib.resources import files as _pkg_files
            cand = str(_pkg_files("fungiforge").joinpath("resources", "af_resistance_panel.tsv"))
            if os.path.exists(cand):
                sys.stderr.write(f"[af_resistance] panel not found at {a.panel}; using packaged copy {cand}\n")
                panel_path = cand
        except Exception:
            pass
    if not os.path.exists(panel_path):
        sys.exit(f"[af_resistance] ERROR: resistance panel not found: {a.panel} (and no packaged copy)")
    curated = read_panel(panel_path)
    fa_rows, evidence, fa_src = ([], {}, None) if a.no_fungamr_panel else fungamr_panel.load_or_derive(a.data_dir, a.fungamr_min_tier)
    panel = merge_panels(curated, fa_rows)
    proteins = read_fasta(a.proteins) if os.path.exists(a.proteins) else {}
    reads = ("sam", read_genotype.load_sam(a.sam)) if a.sam else (a.bam if a.bam and os.path.exists(a.bam) else None)
    genome_depth = read_genotype.genome_mean_depth(a.bam) if (a.bam and os.path.exists(a.bam)) else None
    bundled_ref = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                               "fungiforge", "resources", "af_reference_proteins.faa")
    refidx = load_reference_index(a.data_dir, a.ref_faa, bundled_ref)

    # Illumina base accuracy carries no homopolymer-indel risk, so short-read-only
    # assemblies are high-confidence for point-mutation resistance calls, like hybrid.
    conf = {"hybrid": "high", "illumina_only": "high",
            "ont_only": "provisional_ont_only", "unknown": "provisional"}[a.polish_mode]
    # species relevance: organism_regex names the species ("Aspergillus fumigatus"), the genus
    # ("Candida spp." — the genus match covers it), or "*" / "Fungi" for a pan-fungal row.
    # (A bare "spp." must NOT make a row apply to every genus — caught by the unit tests.)
    genus = species.split()[0] if species and species != "unknown" else ""
    species_resolved = bool(genus) and not re.search(r"\s(sp|spp)\.?$", species)
    def relevant(row):
        org = (row.get("organism_regex", "") or "").strip()
        if not genus:
            return False
        if not row.get("curated", True):          # FungAMR-derived: exact species (or synonym) only
            return species_resolved and re.search(org, species, re.I) is not None
        if org in ("*", "Fungi", "any"):
            return True
        # A genus-level row ("Candida spp.", "Aspergillus") covers every species of the genus. A row
        # that names a species applies to THAT species only: its hotspot numbering and known
        # mutations are that species' (found on DF-005, 2026-09-21: the A. fumigatus cyp51A and
        # hmg1 rows were scanned in an A. flavus isolate through the genus match, and interspecies
        # differences came back as azole resistance). An unresolved "Genus sp." gets genus rows only.
        if re.fullmatch(r"[A-Z][a-z]+(\s+spp?\.?)?", org):
            return re.search(r"\b" + re.escape(genus) + r"\b", org, re.I) is not None
        return species_resolved and re.search(re.escape(species), org, re.I) is not None

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
        # Which reference was used is recorded (`reference`, `reference_match`): a reference from
        # another species carries that species' numbering, so every residue-level result against it
        # is a SCREEN (class cross_species_screen, never counted as resistance), and an unresolved
        # "Genus sp." can only ever be screened. DF-005 (A. flavus, 2026-09-21) was numbered against
        # A. fumigatus cyp51A (77.7 % identity) and hmg1 (85.9 %), and four interspecies differences
        # were reported as azole-resistance calls.
        toks = [t for t in re.split(r"\s+", species.lower()) if len(t) > 2 and t not in ("sp.", "spp.")]
        sp_refs = [(h, s) for (h, _, s) in refs if species_resolved and toks and all(t in h.lower() for t in toks)]
        ref_match = "species"
        if not sp_refs and toks:
            sp_refs = [(h, s) for (h, _, s) in refs if toks[0] in h.lower()]
            ref_match = "genus" if sp_refs else "other"
        if not sp_refs:
            sp_refs = [(h, s) for (h, _, s) in refs]
            ref_match = "other"
        ref_header, refseq = max(sp_refs, key=lambda hs: len(hs[1]))
        cross_species = ref_match != "species"
        orth = best_ortholog(refseq, proteins) if proteins else None
        if not orth:
            calls.append({"gene": gene, "drug_class": row.get("drug_class"), "mechanism": mech,
                          "status": "not_detected", "confidence": conf})
            continue
        name, qseq, pid, aln = orth
        base = {"gene": gene, "drug_class": row.get("drug_class"), "drugs": row.get("drugs"),
                "mechanism": mech, "ortholog": name, "identity": round(pid, 3),
                "confidence": conf, "source": row.get("source", ""),
                "reference": ref_header, "reference_match": ref_match}
        if cross_species:
            base["screen_only"] = True
            base["note"] = (f"reference {ref_header} is not from {species}; residue numbering is species-specific, "
                            "so this is a screen, not a resistance call")
        cds = read_genotype.cds_codon_coords(a.gbk, name, qseq) if (reads is not None and a.gbk and os.path.exists(a.gbk)) else None
        if cds and mech in ("substitution", "overexpression", "GOF") and a.bam and os.path.exists(a.bam):
            lo = min(p for t in cds["codons"] for p in t); hi = max(p for t in cds["codons"] for p in t)
            base["copy_number"] = read_genotype.copy_number(a.bam, cds["contig"], lo, hi, genome_depth)
            if base["copy_number"].get("ratio") and base["copy_number"]["ratio"] >= 1.8:
                base["copy_number"]["flag"] = "possible_duplication"
        if mech == "substitution":
            subs = call_substitutions(gene, refseq, aln, row, evidence, species, curated=row.get("curated", True))
            rs = read_support_for(reads, cds, aln, row, refseq, a.min_reads) if cds else {}
            for c in subs:
                sup = rs.get(c["ref_pos"])
                if sup:
                    sup = dict(sup, agrees=(sup["major"] == c["obs"]) if sup["call"] != "insufficient" else None)
                    c["read_support"] = sup
                c["confidence"] = combined_confidence(a.polish_mode, c.get("read_support"), c.get("evidence_tier"))
                if cross_species:
                    c.update({"class": "cross_species_screen", "known": False, "confidence": "screen_only"})
                    c.pop("note", None)
                calls.append({**base, "status": "variant", **c})
            # known alleles present in the reads but absent from the assembly (minor / heterozygous);
            # meaningless against another species' numbering, so only with a same-species reference
            known, _ = panel_positions(row)
            called_pos = {c["ref_pos"] for c in subs}
            for pos, sup in ({} if cross_species else rs).items():
                if pos in called_pos or sup["call"] == "insufficient":
                    continue
                wt_res = refseq[pos - 1]
                for aa, frac in sup["alleles"].items():
                    if aa == wt_res or aa in ("del", "X", "*"):          # reads agree with the wild-type (or are uninformative)
                        continue
                    if any(mut == aa for (_, mut) in known.get(pos, [])) and frac >= 0.2:
                        change = f"{refseq[pos - 1]}{pos}{aa}"
                        ev = evidence_for(evidence, gene, species, change)
                        ck = row.get("curated_known")
                        tier = "curated" if (row.get("curated", True) and (change in ck if ck is not None else True)) else (ev["tier"] if ev else None)
                        calls.append({**base, "status": "variant", "gene": gene, "change": change, "known": tier != 8,
                                      "class": "known_resistance_mutation" if tier != 8 else "associated_unvalidated",
                                      "evidence_tier": tier, "ref_pos": pos, "wt": refseq[pos - 1], "obs": aa,
                                      "read_support": dict(sup, agrees=None), "assembly_residue": refseq[pos - 1],
                                      "confidence": "high" if frac >= a.min_frac else "provisional_minor_allele",
                                      "note": f"detected in reads only ({sup['call']}, {frac:.0%}); assembly carries the wild-type residue"})
            if not any(c.get("gene") == gene and c.get("status") == "variant" for c in calls):
                wt_call = {**base, "status": "wild_type"}
                if rs:
                    n_ok = sum(1 for v in rs.values() if v["call"] != "insufficient")
                    wt_call["read_support"] = {"positions_checked": len(rs), "positions_covered": n_ok,
                                               "agrees": True if n_ok else None, "min_depth": min(v["depth"] for v in rs.values())}
                    wt_call["confidence"] = combined_confidence(a.polish_mode, wt_call["read_support"])
                calls.append(wt_call)
        elif mech == "loss_of_function":
            trunc = len(qseq) < 0.85 * len(refseq)
            calls.append({**base, "status": "loss_of_function" if trunc else "intact",
                          "note": f"length {len(qseq)}/{len(refseq)} aa" + (" — truncated (candidate LoF)" if trunc else "")})
        elif mech == "GOF":
            subs = call_substitutions(gene, refseq, aln, row)
            calls.append({**base, "status": "present",
                          "gof_variants": [] if cross_species else [c["change"] for c in subs if c["known"]],
                          "note": "regulator; gain-of-function inferred from known residues (confirm with expression/phenotype)"})
        elif mech == "overexpression":
            calls.append({**base, "status": "present",
                          "note": "efflux/target present; over-expression needs RNA-seq — not callable from genome alone"})

    # cyp51A promoter tandem repeats (A. fumigatus)
    tr = None
    if re.search(r"Aspergillus\s+fumigatus", species, re.I) and a.assembly:
        tr = cyp51a_TR.detect_tr(a.assembly, a.gbk)
        if tr.get("tr_site_start") and reads is not None:
            tr["reads"] = tr_read_check(reads, tr, a.min_reads, a.min_frac)
        if tr.get("tr_detected"):
            conf_tr = "high" if a.polish_mode in ("hybrid", "illumina_only") else "medium"
            rd = tr.get("reads") or {}
            if rd.get("call") == "agrees_with_assembly":
                conf_tr = "high"
            elif rd.get("call") in ("reads_lack_copy", "mixed"):
                conf_tr = "discordant"
            calls.append({"gene": "cyp51A_promoter", "drug_class": "azole", "mechanism": "promoter_TR",
                          "status": "resistance", "change": tr["tr_type"], "known": True, "evidence_tier": "curated",
                          "confidence": conf_tr, "read_support": rd or None,
                          "note": tr["note"] + " — pairs with cyp51A L98H (TR34) or Y121F+T289A (TR46)"})
        elif (tr.get("reads") or {}).get("call") == "reads_have_extra_copy":
            rd = tr["reads"]
            L = max(rd["insertion_lengths"], key=rd["insertion_lengths"].get) if rd["insertion_lengths"] else 0
            calls.append({"gene": "cyp51A_promoter", "drug_class": "azole", "mechanism": "promoter_TR",
                          "status": "resistance", "change": "TR34" if 30 <= L <= 40 else "TR46" if 42 <= L <= 52 else f"TR{L}",
                          "known": True, "evidence_tier": "curated", "confidence": "discordant", "read_support": rd,
                          "note": f"reads carry a {L}-bp insertion at the cyp51A TR site that the assembly lacks (tandem repeat missed by assembly)"})

    resistant_classes = sorted({c.get("drug_class") for c in calls
                                if not c.get("screen_only")
                                and (c.get("known") or c.get("status") in ("loss_of_function", "resistance"))
                                and c.get("drug_class")})
    supports = [c["read_support"] for c in calls if isinstance(c.get("read_support"), dict)]
    rs_summary = {"mode": "reads" if reads is not None else "assembly_only",
                  "confirmed": sum(1 for c in calls if c.get("status") == "variant" and (c.get("read_support") or {}).get("agrees") is True),
                  "discordant": sum(1 for c in calls if c.get("confidence") == "discordant"),
                  "insufficient": sum(1 for s_ in supports if s_.get("call") == "insufficient"),
                  "reads_only": sum(1 for c in calls if "reads only" in (c.get("note") or ""))}
    cn_flags = sorted({f"{c['gene']}:{c['copy_number']['ratio']}" for c in calls
                       if isinstance(c.get("copy_number"), dict) and c["copy_number"].get("flag")})
    out = {"sample": a.sample, "stage": "resistance", "species": species,
           "polish_mode": a.polish_mode, "genes_searched": sorted(set(searched)),
           "panel": {"curated_rows": len(curated), "fungamr_rows": len(fa_rows), "fungamr_source": fa_src,
                     "rows_relevant": len([r for r in panel if relevant(r)])},
           "calls": calls, "cyp51A_TR": tr,
           "summary": {"resistant_drug_classes": [c for c in resistant_classes if c],
                       "n_known_mutations": sum(1 for c in calls if c.get("known")),
                       "n_associated_unvalidated": sum(1 for c in calls if c.get("class") == "associated_unvalidated"),
                       "n_associated_other_species": sum(1 for c in calls if c.get("class") == "associated_other_species"),
                       "n_cross_species_screen": sum(1 for c in calls if c.get("class") == "cross_species_screen"),
                       "reference_match": {c["gene"]: c["reference_match"] for c in calls if c.get("reference_match")},
                       "species_resolved": species_resolved,
                       "reference_available": bool(refidx),
                       "read_support": rs_summary, "copy_number_flags": cn_flags}}
    json.dump(out, open(a.out, "w"), indent=2)
    print(f"[af_resistance] {a.sample} · {species} · {len(calls)} calls · "
          f"known={out['summary']['n_known_mutations']} · ref={'yes' if refidx else 'no'}")


if __name__ == "__main__":
    main()
