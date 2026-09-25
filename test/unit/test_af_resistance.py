"""Unit tests for bin/af_resistance.py — panel parsing, ortholog search, substitution / LoF calls,
confidence by polish mode, species relevance, and the cyp51A TR integration (synthetic data)."""
import json, os, subprocess, sys
import pytest

pytest.importorskip("Bio", reason="Biopython is a runtime dependency of the resistance stage (dev extra)")
import af_resistance as afr  # noqa: E402  (sys.path set in conftest)
import cyp51a_TR  # noqa: E402

PANEL = """# synthetic panel
gene\torganism_regex\tdrug_class\tdrugs\tmechanism\thotspot_aa\tknown_mutations\tnote\tsource
cyp51A\tAspergillus fumigatus\tazole\titraconazole,voriconazole\tsubstitution\t54,98,121\tG54W,L98H,Y121F\t-\ttest
cyp51A_promoter\tAspergillus fumigatus\tazole\titraconazole\tpromoter_TR\tNA\tTR34,TR46\t-\ttest
FUR1\tAspergillus spp.\t5-FC\tflucytosine\tloss_of_function\tNA\tNA\t-\ttest
FKS1\tCandida spp.\techinocandin\tcaspofungin\tsubstitution\t639,645\tS639F,S645P\t-\ttest
"""


def make_ref(h):
    """120-aa synthetic cyp51A reference with G at 54, L at 98, Y at 121 (1-based)."""
    p = list(h["random_protein"](130, seed=41))
    p[53], p[97], p[120] = "G", "L", "Y"
    return "".join(p)


def setup(tmp_path, h, ortholog=None, species="Aspergillus fumigatus", extra_prots=None, ref=True, fur1=None):
    panel = tmp_path / "panel.tsv"; panel.write_text(PANEL)
    cyp = make_ref(h)
    fur = h["random_protein"](200, seed=42)
    refs = {"cyp51A__Q4WNT5__Aspergillus_fumigatus": cyp, "FUR1__Q9UVL5__Aspergillus_fumigatus": fur,
            "Cdr1_Erg11_Fcy1_Fks1__P0ABC__Candida_albicans": h["random_protein"](150, seed=43)}
    ref_faa = h["write_fasta"](tmp_path / "ref.faa", refs) if ref else None
    prots = dict(extra_prots or {})
    if ortholog is not None:
        prots["FUN_000123-T1 cyp51A-like"] = ortholog
    if fur1 is not None:
        prots["FUN_000200-T1"] = fur1
    prots["FUN_000999-T1"] = h["random_protein"](300, seed=44)          # unrelated protein
    proteins = h["write_fasta"](tmp_path / "prot.faa", prots)
    sp = tmp_path / "species.txt"; sp.write_text(f"{species}\thigh\n")
    return {"panel": str(panel), "ref": ref_faa, "proteins": proteins, "species": str(sp), "cyp": cyp, "fur": fur}


def run(h, tmp_path, s, polish="hybrid", assembly=None, gbk=None, data_dir="", extra=()):
    out = tmp_path / "res.json"
    cmd = [sys.executable, os.path.join(h["BIN"], "af_resistance.py"), "--sample", "S1", "--proteins", s["proteins"],
           "--species", s["species"], "--panel", s["panel"], "--polish-mode", polish, "--out", str(out), "--data-dir", data_dir, *extra]
    if s["ref"]: cmd += ["--ref-faa", s["ref"]]
    if assembly: cmd += ["--assembly", assembly]
    if gbk: cmd += ["--gbk", gbk]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.load(open(out)), r


def by_gene(doc, gene):
    return [c for c in doc["calls"] if c["gene"] == gene]


def sub(tmp_path, name):
    d = tmp_path / name; d.mkdir(); return d


# ---------------------------------------------------------------- helpers
def test_read_panel_skips_comments_and_maps_columns(tmp_path):
    p = tmp_path / "p.tsv"; p.write_text(PANEL)
    rows = afr.read_panel(p)
    assert [r["gene"] for r in rows] == ["cyp51A", "cyp51A_promoter", "FUR1", "FKS1"]
    assert rows[0]["known_mutations"] == "G54W,L98H,Y121F" and rows[2]["mechanism"] == "loss_of_function"


def test_reference_index_indexes_every_gene_keyword_in_a_header(tmp_path, helpers):
    faa = helpers["write_fasta"](tmp_path / "r.faa", {"Cdr1_Erg11_Fcy1_Fks1__ACC__Candida_albicans": "MKV", "cyp51A__X__Af": "MAA", "CYP51__Y__Ca": "MCC"})
    idx = afr.load_reference_index("", faa, None)
    assert set(idx) >= {"cdr1", "erg11", "fcy1", "fks1", "cyp51a", "cyp51"}
    assert idx["cyp51a"][0][2] == "MAA" and idx["cyp51"][0][2] == "MCC"


def test_best_ortholog_and_residue_mapping(helpers):
    ref = make_ref(helpers)
    q = ref[:97] + "H" + ref[98:]
    prots = {"orth": q, "junk": helpers["random_protein"](130, seed=99)}
    name, seq, pid, aln = afr.best_ortholog(ref, prots)
    assert name == "orth" and pid > 0.99
    assert afr.ref_to_query_residue(aln, 98) == "H" and afr.ref_to_query_residue(aln, 54) == "G"


def test_best_ortholog_returns_none_below_identity_floor(helpers):
    ref = make_ref(helpers)
    assert afr.best_ortholog(ref, {"x": helpers["random_protein"](130, seed=77)}) is None


# ---------------------------------------------------------------- end to end
def test_known_mutation_is_called_with_high_confidence_for_hybrid(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref[:97] + "H" + ref[98:])
    doc, r = run(helpers, tmp_path, s, polish="hybrid")
    c = by_gene(doc, "cyp51A")
    assert len(c) == 1 and c[0]["status"] == "variant" and c[0]["change"] == "L98H" and c[0]["known"] is True
    assert c[0]["class"] == "known_resistance_mutation" and c[0]["confidence"] == "high"
    assert doc["summary"]["resistant_drug_classes"] == ["azole"] and doc["summary"]["n_known_mutations"] == 1
    assert doc["summary"]["reference_available"] is True and "known=1" in r.stdout


def test_wild_type_ortholog_and_novel_hotspot_variant(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref)
    doc, _ = run(helpers, tmp_path, s)
    assert by_gene(doc, "cyp51A")[0]["status"] == "wild_type" and doc["summary"]["resistant_drug_classes"] == []
    b = sub(tmp_path, "b")
    doc2, _ = run(helpers, b, setup(b, helpers, ortholog=ref[:53] + "A" + ref[54:]))
    c = by_gene(doc2, "cyp51A")[0]
    assert c["status"] == "variant" and c["change"] == "G54A" and c["known"] is False and c["class"] == "novel_hotspot_variant"
    assert doc2["summary"]["resistant_drug_classes"] == []                 # novel variants do not claim resistance


def test_ont_only_calls_are_provisional(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref[:97] + "H" + ref[98:])
    doc, _ = run(helpers, tmp_path, s, polish="ont_only")
    assert by_gene(doc, "cyp51A")[0]["confidence"] == "provisional_ont_only" and doc["polish_mode"] == "ont_only"


def test_not_detected_when_no_ortholog(tmp_path, helpers):
    s = setup(tmp_path, helpers, ortholog=None)
    doc, _ = run(helpers, tmp_path, s)
    assert by_gene(doc, "cyp51A")[0]["status"] == "not_detected"


def test_no_reference_is_reported_not_guessed(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref, ref=False)
    doc, _ = run(helpers, tmp_path, s)
    c = by_gene(doc, "cyp51A")[0]
    assert c["status"] == "no_reference" and "install FungAMR" in c["note"]
    assert doc["summary"]["reference_available"] is False


def test_loss_of_function_by_truncation(tmp_path, helpers):
    s = setup(tmp_path, helpers, ortholog=make_ref(helpers))
    fur = s["fur"]
    t = sub(tmp_path, "t")
    doc, _ = run(helpers, t, setup(t, helpers, ortholog=make_ref(helpers), fur1=fur[:90]))
    c = by_gene(doc, "FUR1")[0]
    assert c["status"] == "loss_of_function" and "truncated" in c["note"]
    assert "5-FC" in doc["summary"]["resistant_drug_classes"]
    f = sub(tmp_path, "f")
    doc2, _ = run(helpers, f, setup(f, helpers, ortholog=make_ref(helpers), fur1=fur))
    assert by_gene(doc2, "FUR1")[0]["status"] == "intact"


def test_species_relevance_filters_panel_rows(tmp_path, helpers):
    s = setup(tmp_path, helpers, ortholog=make_ref(helpers))
    doc, _ = run(helpers, tmp_path, s)
    assert doc["genes_searched"] == ["FUR1", "cyp51A", "cyp51A_promoter"]      # Candida FKS1 row not searched
    c = sub(tmp_path, "c")
    doc2, _ = run(helpers, c, setup(c, helpers, species="Candida auris"))
    assert doc2["genes_searched"] == ["FKS1"] and doc2["cyp51A_TR"] is None
    u = sub(tmp_path, "u")
    doc3, _ = run(helpers, u, setup(u, helpers, species="unknown"))
    assert doc3["genes_searched"] == [] and doc3["calls"] == []


def test_tr34_detected_via_assembly_and_gbk(tmp_path, helpers):
    unit = helpers["TR34_UNIT"]
    promoter = helpers["random_dna"](532, seed=51) + unit * 2
    gene = helpers["random_dna"](900, seed=52)
    seq = helpers["random_dna"](100, seed=53) + promoter + gene + helpers["random_dna"](100, seed=54)
    gs = 100 + len(promoter)
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(gs, gs + len(gene), 1, "gene", {"gene": ["cyp51A"]})])
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    s = setup(tmp_path, helpers, ortholog=make_ref(helpers))
    doc, _ = run(helpers, tmp_path, s, polish="hybrid", assembly=asm, gbk=gbk)
    tr = by_gene(doc, "cyp51A_promoter")
    assert len(tr) == 1 and tr[0]["status"] == "resistance" and tr[0]["change"] == "TR34" and tr[0]["confidence"] == "high"
    assert doc["cyp51A_TR"]["tr_type"] == "TR34" and doc["summary"]["resistant_drug_classes"] == ["azole"]
    doc_ont, _ = run(helpers, tmp_path, s, polish="ont_only", assembly=asm, gbk=gbk)
    assert by_gene(doc_ont, "cyp51A_promoter")[0]["confidence"] == "medium"


def test_missing_panel_fails_loudly(tmp_path, helpers):
    s = setup(tmp_path, helpers, ortholog=make_ref(helpers))
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "af_resistance.py"), "--sample", "S1", "--proteins", s["proteins"],
                        "--species", s["species"], "--panel", str(tmp_path / "missing.tsv"), "--out", str(tmp_path / "o.json")],
                       capture_output=True, text=True)
    # the packaged panel is a legitimate fallback when the fungiforge package is importable; otherwise it must exit non-zero
    assert r.returncode != 0 or "using packaged copy" in r.stderr


# ---------------------------------------------------------------- W2.4: FungAMR evidence, read support, TR reads
import read_genotype as rg  # noqa: E402
from test_read_genotype import encode, sam_read  # noqa: E402

FUNGAMR_HDR = ["species", "gene or protein name", "mutation", "drug", "confidence score", "pubmedid"]


def fungamr_dir(tmp_path, *rows):
    d = tmp_path / "db" / "fungamr"; d.mkdir(parents=True)
    (d / "FungAMR_070425.tsv").write_text("\t".join(FUNGAMR_HDR) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
    return str(tmp_path / "db")


def test_fungamr_rows_merge_into_the_curated_row_and_add_tiers(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref[:53] + "W" + ref[54:97] + "H" + ref[98:])       # G54W (curated) + L98H
    alt46 = "Y" if ref[45] != "Y" else "F"
    m46 = f"{ref[45]}46{alt46}"
    db = fungamr_dir(tmp_path, ["Aspergillus fumigatus", "cyp51A", "L98H", "Itraconazole", "1", "1"],
                     ["Aspergillus fumigatus", "cyp51A", "G54W", "Itraconazole", "8", "2"],
                     ["Aspergillus fumigatus", "cyp51A", m46, "Itraconazole", "8", "3"])     # tier-8 only, not in the curated row
    (tmp_path / "b").mkdir()
    s2 = setup(tmp_path / "b", helpers, ortholog=ref[:45] + alt46 + ref[46:])
    doc, _ = run(helpers, tmp_path, s, data_dir=db)
    calls = {c["change"]: c for c in by_gene(doc, "cyp51A")}
    assert calls["L98H"]["evidence_tier"] == "curated" and calls["L98H"]["known"] and calls["G54W"]["evidence_tier"] == "curated"
    assert doc["panel"]["fungamr_rows"] == 1 and doc["panel"]["fungamr_source"] == "FungAMR_070425.tsv"
    # the FungAMR-only F46Y is reported as an association (tier 8), not as resistance
    doc2, _ = run(helpers, tmp_path / "b", s2, data_dir=db)
    c = by_gene(doc2, "cyp51A")
    assert len(c) == 1 and c[0]["change"] == m46 and c[0]["class"] == "associated_unvalidated" and c[0]["known"] is False
    assert c[0]["evidence_tier"] == 8 and c[0]["confidence"] == "low_evidence"
    assert doc2["summary"]["resistant_drug_classes"] == [] and doc2["summary"]["n_associated_unvalidated"] == 1
    doc3, _ = run(helpers, tmp_path / "b", s2, data_dir=db, extra=["--no-fungamr-panel"])
    assert by_gene(doc3, "cyp51A")[0]["status"] == "wild_type"


def test_merge_panels_and_evidence_lookup():
    cur = [{"gene": "ERG11", "organism_regex": "Candida albicans", "mechanism": "substitution", "hotspot_aa": "132", "known_mutations": "Y132F"}]
    fa = [{"gene": "Erg11", "organism_regex": r"Candida\ albicans", "mechanism": "substitution", "hotspot_aa": "132,143", "known_mutations": "Y132F,K143R", "source": "FungAMR 1"},
          {"gene": "Erg11", "organism_regex": r"Candidozyma\ auris|Candida\ auris", "mechanism": "substitution", "hotspot_aa": "132", "known_mutations": "Y132F", "source": "FungAMR 2"}]
    m = afr.merge_panels(cur, fa)
    assert len(m) == 2 and m[0]["known_mutations"] == "Y132F,K143R" and m[0]["hotspot_aa"] == "132,143" and m[0]["fungamr_merged"] and m[0]["curated"]
    assert m[1]["curated"] is False and m[1]["gene"] == "Erg11"
    ev = {("erg11", "candida auris", "Y132F"): {"tier": 2, "classes": ["azole"]}}
    assert afr.evidence_for(ev, "ERG11", "Candida auris", "Y132F")["tier"] == 2
    assert afr.evidence_for(ev, "ERG11", "Candida haemulonii", "Y132F")["tier"] == 2      # genus fallback
    assert afr.evidence_for(ev, "FKS1", "Candida auris", "Y132F") is None


def test_combined_confidence_rules():
    assert afr.combined_confidence("hybrid") == "high" and afr.combined_confidence("ont_only") == "provisional_ont_only"
    assert afr.combined_confidence("ont_only", {"agrees": True}) == "high"
    assert afr.combined_confidence("hybrid", {"agrees": False}) == "discordant"
    assert afr.combined_confidence("hybrid", {"agrees": None}) == "high"
    assert afr.combined_confidence("hybrid", {"agrees": True}, tier=8) == "low_evidence"


def gene_locus(helpers, protein, promoter_units=1):
    """Contig with a wild-type-style cyp51A promoter (one or two TR34 units) and a single-exon
    CDS encoding `protein`; returns (seq, cds_start0, cds_end0, tr_unit_offset0)."""
    unit = helpers["TR34_UNIT"]
    left = helpers["random_dna"](400, seed=61)
    promoter = left + unit * promoter_units + helpers["random_dna"](250, seed=62)
    cds = encode(protein) + "TAA"
    seq = helpers["random_dna"](100, seed=63) + promoter + cds + helpers["random_dna"](300, seed=64)
    cs = 100 + len(promoter)
    return seq, cs, cs + len(cds), 100 + len(left)


def reads_from(seq, start, end, n=40, length=250, variants=None, ins=None, dele=None):
    """Reads tiling seq[start:end] (0-based); `variants` = {pos0: base} applied to every even read;
    `ins` = (pos0, seq) insertion / `dele` = (pos0, len) deletion in every read."""
    out = []
    span = max(1, end - start - length)
    for i in range(n):
        p0 = start + (i * 13) % span
        s = list(seq[p0:p0 + length])
        if variants and i % 2 == 0:
            for pos0, b in variants.items():
                if p0 <= pos0 < p0 + length:
                    s[pos0 - p0] = b
        s = "".join(s)
        cigar = f"{len(s)}M"
        if ins and p0 + 30 <= ins[0] < p0 + length - 30:
            k = ins[0] - p0
            s, cigar = s[:k] + ins[1] + s[k:], f"{k}M{len(ins[1])}I{length - k}M"
        if dele and p0 + 30 <= dele[0] < p0 + length - 30 - dele[1]:
            k = dele[0] - p0
            s, cigar = s[:k] + s[k + dele[1]:], f"{k}M{dele[1]}D{length - k - dele[1]}M"
        out.append(sam_read(f"r{i}", p0 + 1, s, cigar=cigar))
    return out


def run_reads(h, tmp_path, s, sam_lines, assembly, gbk, polish="ont_only"):
    sam = tmp_path / "reads.sam"; sam.write_text("@HD\tVN:1.6\n" + "\n".join(sam_lines) + "\n")
    return run(h, tmp_path, s, polish=polish, assembly=assembly, gbk=gbk, extra=["--sam", str(sam)])


def test_read_support_confirms_an_assembly_variant_and_raises_ont_confidence(tmp_path, helpers):
    ref = make_ref(helpers)
    mut = ref[:97] + "H" + ref[98:]
    seq, cs, ce, _ = gene_locus(helpers, mut)
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(cs, ce, 1, "gene", {"gene": ["cyp51A"]}),
                                                         (cs, ce, 1, "CDS", {"locus_tag": ["FUN_000123"], "protein_id": ["ncbi:FUN_000123-T1"], "translation": [mut]})])
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    s = setup(tmp_path, helpers, ortholog=mut)
    doc, _ = run_reads(helpers, tmp_path, s, reads_from(seq, cs - 400, ce + 200, n=80), asm, gbk, polish="ont_only")   # reads also span the TR site
    c = by_gene(doc, "cyp51A")[0]
    assert c["change"] == "L98H" and c["read_support"]["call"] == "homozygous" and c["read_support"]["agrees"] is True
    assert c["read_support"]["major"] == "H" and c["read_support"]["depth"] >= 10 and c["confidence"] == "high"
    assert doc["summary"]["read_support"] == {"mode": "reads", "confirmed": 1, "discordant": 0, "insufficient": 0, "reads_only": 0}
    # the TR site of the wild-type promoter: reads agree (no indel), no TR call
    assert doc["cyp51A_TR"]["reads"]["call"] == "agrees_with_assembly" and not doc["cyp51A_TR"]["tr_detected"]
    # reads that contradict the assembly (all wild-type) make the call discordant
    doc2, _ = run_reads(helpers, tmp_path, s, reads_from(seq.replace(encode(mut), encode(ref)), cs - 400, ce + 200, n=80), asm, gbk, polish="hybrid")
    assert by_gene(doc2, "cyp51A")[0]["confidence"] == "discordant"


def test_heterozygous_known_allele_in_reads_only_is_reported(tmp_path, helpers):
    ref = make_ref(helpers)
    seq, cs, ce, _ = gene_locus(helpers, ref)
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(cs, ce, 1, "CDS", {"locus_tag": ["FUN_000123"], "translation": [ref]})])
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    s = setup(tmp_path, helpers, ortholog=ref)
    codon98 = cs + 97 * 3
    hcod = encode("H")
    variants = {codon98 + k: hcod[k] for k in range(3)}
    doc, _ = run_reads(helpers, tmp_path, s, reads_from(seq, cs - 100, ce + 100, n=60, variants=variants), asm, gbk, polish="hybrid")
    c = [x for x in by_gene(doc, "cyp51A") if x["status"] == "variant"]
    assert len(c) == 1 and c[0]["change"] == "L98H" and "reads only" in c[0]["note"] and c[0]["read_support"]["call"] == "heterozygous"
    assert c[0]["confidence"] == "provisional_minor_allele" and c[0]["assembly_residue"] == "L" and doc["summary"]["read_support"]["reads_only"] == 1
    assert doc["summary"]["resistant_drug_classes"] == ["azole"]
    # a known-mutation token whose mutant letter equals the reference residue (FungAMR numbering
    # drift) must not turn wild-type reads into a "reads only" call (found in the image check)
    db = fungamr_dir(tmp_path, ["Aspergillus fumigatus", "cyp51A", f"K54{ref[53]}", "Itraconazole", "1", "1"])
    doc2, _ = run_reads(helpers, tmp_path, s, reads_from(seq, cs - 100, ce + 100, n=60), asm, gbk, polish="hybrid")
    assert not [x for x in by_gene(doc2, "cyp51A") if x["status"] == "variant"]
    docdb, _ = run(helpers, tmp_path, s, polish="hybrid", assembly=asm, gbk=gbk, data_dir=db, extra=["--sam", str(tmp_path / "reads.sam")])
    assert not [x for x in by_gene(docdb, "cyp51A") if x["status"] == "variant"] and by_gene(docdb, "cyp51A")[0]["read_support"]["agrees"] is True


def test_tr34_in_reads_but_not_in_assembly_and_the_reverse(tmp_path, helpers):
    ref = make_ref(helpers)
    unit = helpers["TR34_UNIT"]
    seq, cs, ce, u0 = gene_locus(helpers, ref, promoter_units=1)
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(cs, ce, 1, "gene", {"gene": ["cyp51A"]}), (cs, ce, 1, "CDS", {"locus_tag": ["FUN_000123"], "translation": [ref]})])
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    s = setup(tmp_path, helpers, ortholog=ref)
    doc, _ = run_reads(helpers, tmp_path, s, reads_from(seq, u0 - 300, u0 + 400, n=40, ins=(u0 + 34, unit)), asm, gbk, polish="hybrid")
    tr = doc["cyp51A_TR"]
    assert tr["tr_detected"] is False and tr["tr_site_start"] == u0 + 1 and tr["reads"]["call"] == "reads_have_extra_copy" and tr["reads"]["insertion_lengths"] == {"34": tr["reads"]["with_insertion"]}
    c = by_gene(doc, "cyp51A_promoter")
    assert len(c) == 1 and c[0]["change"] == "TR34" and c[0]["confidence"] == "discordant" and "assembly lacks" in c[0]["note"]
    assert doc["summary"]["resistant_drug_classes"] == ["azole"]
    # assembly with two copies, reads lacking one -> discordant; reads agreeing -> high
    seq2, cs2, ce2, u2 = gene_locus(helpers, ref, promoter_units=2)
    gbk2 = helpers["write_gbk"](tmp_path / "g2.gbk", seq2, [(cs2, ce2, 1, "gene", {"gene": ["cyp51A"]}), (cs2, ce2, 1, "CDS", {"locus_tag": ["FUN_000123"], "translation": [ref]})])
    asm2 = helpers["write_fasta"](tmp_path / "a2.fa", {"contig_1": seq2})
    d_ok, _ = run_reads(helpers, tmp_path, s, reads_from(seq2, u2 - 300, u2 + 500, n=40), asm2, gbk2, polish="ont_only")
    assert d_ok["cyp51A_TR"]["tr_type"] == "TR34" and d_ok["cyp51A_TR"]["reads"]["call"] == "agrees_with_assembly" and by_gene(d_ok, "cyp51A_promoter")[0]["confidence"] == "high"
    d_bad, _ = run_reads(helpers, tmp_path, s, reads_from(seq2, u2 - 300, u2 + 500, n=40, dele=(u2 + 34, 34)), asm2, gbk2, polish="hybrid")
    assert d_bad["cyp51A_TR"]["reads"]["call"] == "reads_lack_copy" and by_gene(d_bad, "cyp51A_promoter")[0]["confidence"] == "discordant"


# ---------------------------------------------------------------- DF-005 regression (2026-09-25)
GENUS_PANEL = """# genus-level row + species-named rows
gene\torganism_regex\tdrug_class\tdrugs\tmechanism\thotspot_aa\tknown_mutations\tnote\tsource
cyp51A\tAspergillus spp.\tazole\titraconazole\tsubstitution\t54,98,121\tG54W,L98H,Y121F\t-\ttest
hmg1\tAspergillus fumigatus\tazole\titraconazole\tsubstitution\t995\tV995I\t-\ttest
"""


def test_species_named_row_applies_to_that_species_only(tmp_path, helpers):
    """DF-005: the A. fumigatus cyp51A/hmg1 rows were scanned in an A. flavus isolate through the genus
    match and interspecies differences became azole-resistance calls. A species-named row is for
    that species; a genus row ('Aspergillus spp.') covers the genus; 'Aspergillus sp.' gets genus rows only."""
    ref = make_ref(helpers)
    for species, expected in (("Aspergillus flavus", ["cyp51A"]), ("Aspergillus sp.", ["cyp51A"]),
                              ("Aspergillus fumigatus", ["cyp51A", "hmg1"])):
        d = sub(tmp_path, species.replace(" ", "_").replace(".", ""))
        s = setup(d, helpers, ortholog=ref, species=species)
        (d / "panel.tsv").write_text(GENUS_PANEL)
        doc, _ = run(helpers, d, s)
        assert doc["genes_searched"] == expected, (species, doc["genes_searched"])


def test_other_species_reference_makes_every_residue_call_a_screen(tmp_path, helpers):
    """Only an A. fumigatus cyp51A reference is staged and the isolate is A. flavus (genus row applies):
    the L98H-shaped difference is reported, but as a cross-species screen — not known, not counted."""
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref[:97] + "H" + ref[98:], species="Aspergillus flavus")
    (tmp_path / "panel.tsv").write_text(GENUS_PANEL)
    doc, r = run(helpers, tmp_path, s)
    c = by_gene(doc, "cyp51A")
    assert len(c) == 1 and c[0]["change"] == "L98H" and c[0]["known"] is False
    assert c[0]["class"] == "cross_species_screen" and c[0]["confidence"] == "screen_only" and c[0]["screen_only"] is True
    assert c[0]["reference"] == "cyp51A__Q4WNT5__Aspergillus_fumigatus" and c[0]["reference_match"] == "genus"
    assert "not from Aspergillus flavus" in c[0]["note"]
    assert doc["summary"]["resistant_drug_classes"] == [] and doc["summary"]["n_known_mutations"] == 0
    assert doc["summary"]["n_cross_species_screen"] == 1 and doc["summary"]["reference_match"] == {"cyp51A": "genus"}
    assert doc["summary"]["species_resolved"] is True and "known=0" in r.stdout


def test_same_species_reference_is_a_real_call_and_is_recorded(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref[:97] + "H" + ref[98:])
    doc, _ = run(helpers, tmp_path, s)
    c = by_gene(doc, "cyp51A")[0]
    assert c["known"] is True and c["reference_match"] == "species" and "screen_only" not in c
    assert c["reference"] == "cyp51A__Q4WNT5__Aspergillus_fumigatus" and doc["summary"]["resistant_drug_classes"] == ["azole"]


def test_unresolved_genus_only_species_is_screened_never_called(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref[:97] + "H" + ref[98:], species="Aspergillus sp.")
    (tmp_path / "panel.tsv").write_text(GENUS_PANEL)
    doc, _ = run(helpers, tmp_path, s)
    c = by_gene(doc, "cyp51A")
    assert c and all(x["class"] == "cross_species_screen" for x in c) and doc["summary"]["species_resolved"] is False
    assert doc["summary"]["resistant_drug_classes"] == [] and doc["summary"]["reference_match"]["cyp51A"] == "genus"


def test_evidence_from_another_species_of_the_genus_is_an_association_only():
    ev = {("hmg1", "aspergillus fumigatus", "V995I"): {"tier": 1, "classes": ["azole"]}}
    row = {"gene": "hmg1", "hotspot_aa": "995", "known_mutations": "V995I", "curated": False}
    ref = "MKLVSTQEWRYPGA" * 71 + "V" + "MKLVS"          # 1000 aa, V at 995
    qry = ref[:994] + "I" + ref[995:]
    calls = afr.call_substitutions("hmg1", ref, _aln(ref, qry), row, evidence=ev, species="Aspergillus flavus", curated=False)
    assert len(calls) == 1 and calls[0]["class"] == "associated_other_species" and calls[0]["known"] is False
    assert calls[0]["evidence_species"] == "aspergillus fumigatus"
    same = afr.call_substitutions("hmg1", ref, _aln(ref, qry), row, evidence=ev, species="Aspergillus fumigatus", curated=False)
    assert same[0]["class"] == "known_resistance_mutation" and same[0]["known"] is True


def _aln(ref, qry):
    """The alignment object call_substitutions expects, built the way best_ortholog builds it."""
    return afr.best_ortholog(ref, {"q": qry})[3]
