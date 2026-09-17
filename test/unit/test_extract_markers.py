"""Unit tests for bin/extract_markers.py — tblastn HSP chaining across introns, strand handling,
LSU D1/D2 from barrnap coordinates, and the CLI contract (one FASTA per locus + JSON summary)."""
import json, os, subprocess, sys
import extract_markers as em  # noqa: E402  (sys.path set in conftest)


def tbl(tmp_path, rows):
    """rows: (qseqid, sseqid, pident, length, qstart, qend, sstart, send, evalue, bitscore)"""
    p = tmp_path / "tblastn.tsv"
    p.write_text("".join("\t".join(str(x) for x in r) + "\n" for r in rows))
    return str(p)


def test_parse_tblastn_groups_by_marker_and_orients_strand(tmp_path):
    p = tbl(tmp_path, [("CaM|XP_1|Af", "c1", 95.0, 100, 1, 100, 1000, 1299, 1e-50, 200.0),
                       ("BenA|XP_2|Af", "c2", 90.0, 80, 1, 80, 5240, 5001, 1e-40, 150.0), ("junk",)])
    h = em.parse_tblastn(p)
    assert set(h) == {"CaM", "BenA"}
    assert h["CaM"][0]["strand"] == "+" and h["CaM"][0]["sstart"] == 1000 and h["CaM"][0]["send"] == 1299
    assert h["BenA"][0]["strand"] == "-" and h["BenA"][0]["sstart"] == 5001 and h["BenA"][0]["send"] == 5240
    assert em.parse_tblastn(str(tmp_path / "none.tsv")) == {} and em.parse_tblastn(None) == {}


def test_chain_best_joins_exons_within_gap_and_picks_best_contig():
    hs = [{"contig": "c1", "strand": "+", "sstart": 1000, "send": 1300, "qstart": 1, "qend": 100, "pident": 95, "bitscore": 200},
          {"contig": "c1", "strand": "+", "sstart": 1400, "send": 1700, "qstart": 101, "qend": 200, "pident": 90, "bitscore": 180},   # exon 2 after a 100-bp intron
          {"contig": "c1", "strand": "+", "sstart": 90000, "send": 90300, "qstart": 1, "qend": 100, "pident": 60, "bitscore": 90},    # paralog far away
          {"contig": "c7", "strand": "-", "sstart": 10, "send": 310, "qstart": 1, "qend": 100, "pident": 70, "bitscore": 120}]
    loc = em.chain_best(hs, max_gap=5000, min_bitscore=80)
    assert (loc["contig"], loc["strand"], loc["start"], loc["end"], loc["n_hsps"]) == ("c1", "+", 1000, 1700, 2)
    assert loc["bitscore"] == 380 and loc["query_residues_covered"] == 200 and loc["n_contigs_hit"] == 2
    assert em.chain_best(hs, min_bitscore=500) is None and em.chain_best([]) is None
    # the far paralog is not chained; with a huge gap it would be
    assert em.chain_best(hs, max_gap=100000, min_bitscore=80)["n_hsps"] == 3


def test_lsu_from_gff_takes_longest_28s_strand_aware(tmp_path, helpers):
    seq = helpers["random_dna"](5000, seed=7)
    contigs = {"c1": seq, "c2": "A" * 100}
    gff = tmp_path / "rrna.gff"
    gff.write_text("##gff-version 3\n"
                   "c1\tbarrnap\trRNA\t100\t400\t0\t+\t.\tName=18S_rRNA;product=18S ribosomal RNA (partial)\n"
                   "c1\tbarrnap\trRNA\t1000\t4400\t0\t-\t.\tName=28S_rRNA;product=28S ribosomal RNA\n"
                   "c1\tbarrnap\trRNA\t4500\t4700\t0\t+\t.\tName=28S_rRNA;product=28S ribosomal RNA (partial)\n"
                   "cX\tbarrnap\trRNA\t1\t9000\t0\t+\t.\tName=28S_rRNA;product=28S ribosomal RNA\n")
    lsu = em.lsu_from_gff(str(gff), contigs, lsu_len=900)
    assert (lsu["contig"], lsu["start"], lsu["end"], lsu["strand"]) == ("c1", 1000, 4400, "-")
    assert lsu["seq"] == em.revcomp(seq[999:4400])[:900] and lsu["full_len"] == 3401 and lsu["truncated"] is False
    assert em.lsu_from_gff(str(tmp_path / "none.gff"), contigs) is None and em.lsu_from_gff(None, contigs) is None
    gff.write_text("c1\tbarrnap\trRNA\t1000\t1300\t0\t+\t.\tName=28S_rRNA;product=28S ribosomal RNA (partial)\n")
    short = em.lsu_from_gff(str(gff), contigs, lsu_len=900)
    assert short["truncated"] is True and len(short["seq"]) == 301 and short["partial"] is True


def run(tmp_path, h, fasta, tblastn=None, gff=None, extra=()):
    cmd = [sys.executable, os.path.join(h["BIN"], "extract_markers.py"), "--sample", "S1", "--assembly", fasta,
           "--out-prefix", str(tmp_path / "S1.marker"), "--json", str(tmp_path / "m.json"), *extra]
    if tblastn: cmd += ["--tblastn", tblastn]
    if gff: cmd += ["--rrna-gff", gff]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.load(open(tmp_path / "m.json")), r.stdout


def read_fa(p):
    lines = open(p).read().splitlines()
    return (lines[0], "".join(lines[1:])) if lines else (None, "")


def test_cli_extracts_each_locus_with_flank_in_coding_orientation(tmp_path, helpers):
    seq = helpers["random_dna"](6000, seed=3)
    fa = helpers["write_fasta"](tmp_path / "asm.fa", {"c1": seq}, width=80)
    t = tbl(tmp_path, [("CaM|XP_1|Af", "c1", 95.0, 100, 1, 100, 1001, 1300, 1e-50, 200.0),
                       ("CaM|XP_1|Af", "c1", 90.0, 100, 101, 200, 1401, 1700, 1e-40, 180.0),
                       ("BenA|XP_2|Af", "c1", 92.0, 100, 1, 100, 4300, 4001, 1e-45, 190.0)])
    gff = tmp_path / "rrna.gff"
    gff.write_text("c1\tbarrnap\trRNA\t5000\t5950\t0\t+\t.\tName=28S_rRNA;product=28S ribosomal RNA\n")
    doc, out = run(tmp_path, helpers, fa, t, str(gff), extra=["--flank", "50"])
    cam = doc["loci"]["CaM"]
    assert cam["found"] and (cam["extracted_start"], cam["extracted_end"], cam["strand"], cam["n_hsps"]) == (951, 1750, "+", 2)
    hdr, s = read_fa(tmp_path / "S1.marker.CaM.fasta")
    assert hdr.startswith(">S1|CaM|c1:951-1750(+)") and s == seq[950:1750] and cam["length"] == 800
    hdr2, s2 = read_fa(tmp_path / "S1.marker.BenA.fasta")
    assert "(-)" in hdr2 and s2 == em.revcomp(seq[3950:4350])                     # minus strand -> reverse-complemented
    assert doc["loci"]["TEF1"] == {"found": False} and os.path.getsize(tmp_path / "S1.marker.TEF1.fasta") == 0
    assert doc["loci"]["RPB2"]["found"] is False
    lsu = doc["loci"]["LSU"]
    assert lsu["found"] and lsu["length"] == 900 and lsu["full_28S_len"] == 951 and lsu["truncated"] is False
    assert read_fa(tmp_path / "S1.marker.LSU.fasta")[1] == seq[4999:5899]
    assert "found CaM, BenA, LSU" in out and "missing TEF1, RPB2" in out


def test_cli_without_hits_writes_empty_files_and_summary(tmp_path, helpers):
    fa = helpers["write_fasta"](tmp_path / "asm.fa", {"c1": "ACGT" * 100})
    doc, out = run(tmp_path, helpers, fa)
    assert all(v == {"found": False} for v in doc["loci"].values()) and set(doc["loci"]) == {"CaM", "BenA", "TEF1", "RPB2", "LSU"}
    assert all(os.path.getsize(tmp_path / f"S1.marker.{m}.fasta") == 0 for m in doc["loci"])
    assert "found none" in out


def test_flank_is_clipped_at_contig_ends(tmp_path, helpers):
    fa = helpers["write_fasta"](tmp_path / "asm.fa", {"c1": helpers["random_dna"](500, seed=1)})
    t = tbl(tmp_path, [("TEF1|XP_3|Af", "c1", 95.0, 100, 1, 100, 10, 309, 1e-50, 200.0)])
    doc, _ = run(tmp_path, helpers, fa, t)
    assert (doc["loci"]["TEF1"]["extracted_start"], doc["loci"]["TEF1"]["extracted_end"]) == (1, 500)
