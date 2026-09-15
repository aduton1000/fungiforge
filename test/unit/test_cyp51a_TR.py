"""Unit tests for bin/cyp51a_TR.py — the structural cyp51A promoter tandem-repeat detector."""
import json, os, subprocess, sys
import pytest

pytest.importorskip("Bio", reason="Biopython is a runtime dependency of the resistance stage (dev extra)")
import cyp51a_TR as tr  # noqa: E402  (sys.path set in conftest)


def genome_with_promoter(h, tr_copies, gene_len=1500, upstream=600, seed=3):
    """Return (seq, gene_start0, gene_end0): a plus-strand gene whose 600-bp upstream promoter
    ends with `tr_copies` tandem copies of the TR34 unit (0 copies = wild-type promoter)."""
    unit = h["TR34_UNIT"]
    filler = h["random_dna"](upstream - tr_copies * len(unit), seed=seed)
    promoter = filler + unit * tr_copies
    gene = h["random_dna"](gene_len, seed=seed + 10)
    tail = h["random_dna"](300, seed=seed + 20)
    seq = tail + promoter + gene + tail
    start = len(tail) + len(promoter)
    return seq, start, start + gene_len


def test_tandem_scan_finds_two_copies_of_a_34bp_unit(helpers):
    promoter = helpers["random_dna"](400, seed=5) + helpers["TR34_UNIT"] * 2 + helpers["random_dna"](50, seed=6)
    L, copies, pos = tr._tandem_scan(promoter)
    # the scanner may start a few bases early when the filler happens to end like the unit (a rotation)
    assert (L, copies) == (34, 2) and 390 <= pos <= 400


def test_tandem_scan_three_copies_and_tolerates_mismatch(helpers):
    unit = helpers["TR34_UNIT"]
    mutated = unit[:5] + ("A" if unit[5] != "A" else "C") + unit[6:]   # 33/34 identical (97 %)
    promoter = helpers["random_dna"](200, seed=7) + unit + mutated + unit
    L, copies, pos = tr._tandem_scan(promoter)
    assert L == 34 and copies == 3 and 190 <= pos <= 200


def test_tandem_scan_wild_type_promoter_is_none(helpers):
    assert tr._tandem_scan(helpers["random_dna"](600, seed=11)) is None


def test_revcomp():
    assert tr._revcomp("ACGTN") == "NACGT" and tr._revcomp("AAC") == "GTT"


def test_detect_tr_plus_strand_named_gene(tmp_path, helpers):
    seq, s, e = genome_with_promoter(helpers, 2)
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(s, e, 1, "gene", {"gene": ["cyp51A"], "locus_tag": ["FUN_001"]})])
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    res = tr.detect_tr(asm, gbk)
    assert res["cyp51A_located"] and res["tr_detected"] and res["tr_type"] == "TR34"
    assert res["unit_len"] == 34 and res["copies"] == 2 and "2×34 bp" in res["note"]


def test_detect_tr_wild_type_promoter(tmp_path, helpers):
    seq, s, e = genome_with_promoter(helpers, 0)
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(s, e, 1, "gene", {"gene": ["cyp51A"]})])
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    res = tr.detect_tr(asm, gbk)
    assert res["cyp51A_located"] and not res["tr_detected"] and res["tr_type"] == "none"
    assert "wild-type promoter" in res["note"]


def test_detect_tr_minus_strand_gene(tmp_path, helpers):
    unit = helpers["TR34_UNIT"]
    promoter = helpers["random_dna"](532, seed=21) + unit * 2               # 600 bp, ends at the gene start
    gene = helpers["random_dna"](1200, seed=22)
    plus_layout = promoter + gene                                             # what the gene "sees" on its own strand
    seq = helpers["random_dna"](200, seed=23) + tr._revcomp(plus_layout) + helpers["random_dna"](200, seed=24)
    gs, ge = 200, 200 + len(gene)                                             # gene occupies the first part of the revcomp block
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(gs, ge, -1, "CDS", {"gene": ["erg11"], "product": ["cyp51A"]})])
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    res = tr.detect_tr(asm, gbk)
    assert res["tr_detected"] and res["tr_type"] == "TR34" and res["copies"] == 2


def test_detect_tr_prefers_named_cyp51a_over_cyp51b_paralog_by_product(tmp_path, helpers):
    seq, s, e = genome_with_promoter(helpers, 2)
    # a CYP51B paralog with the demethylase product wording appears FIRST, then the real cyp51A
    feats = [(50, 250, 1, "CDS", {"product": ["Sterol 14-alpha demethylase CYP51B"]}),
             (s, e, 1, "gene", {"gene": ["cyp51A"]})]
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, feats)
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    assert tr.detect_tr(asm, gbk)["tr_detected"]


def test_detect_tr_product_fallback_when_no_gene_name(tmp_path, helpers):
    seq, s, e = genome_with_promoter(helpers, 2)
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(s, e, 1, "CDS", {"product": ["eburicol 14-alpha-demethylase"]})])
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    assert tr.detect_tr(asm, gbk)["tr_type"] == "TR34"


def test_detect_tr_without_gbk_or_without_cyp51a(tmp_path, helpers):
    seq, s, e = genome_with_promoter(helpers, 2)
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    res = tr.detect_tr(asm, None)
    assert not res["cyp51A_located"] and "TR scan skipped" in res["note"]
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(s, e, 1, "gene", {"gene": ["fks1"]})])
    assert not tr.detect_tr(asm, gbk)["cyp51A_located"]


def test_tr46_triplicate_naming(helpers):
    unit46 = helpers["random_dna"](46, seed=31)
    promoter = helpers["random_dna"](300, seed=32) + unit46 * 3
    L, copies, _ = tr._tandem_scan(promoter)
    assert (L, copies) == (46, 3)


def test_cli_writes_json(tmp_path, helpers):
    seq, s, e = genome_with_promoter(helpers, 2)
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(s, e, 1, "gene", {"gene": ["cyp51A"]})])
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "cyp51a_TR.py"), "--assembly", asm, "--gbk", gbk,
                        "--out", str(tmp_path / "tr.json")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert json.load(open(tmp_path / "tr.json"))["tr_type"] == "TR34" and '"tr_detected": true' in r.stdout
