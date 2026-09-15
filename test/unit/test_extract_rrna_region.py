"""Unit tests for bin/extract_rrna_region.py — padded, merged rRNA-operon windows from barrnap GFF."""
import os, subprocess, sys
import extract_rrna_region as er  # noqa: E402  (sys.path set in conftest)


def gff(tmp_path, rows):
    p = tmp_path / "rrna.gff"
    p.write_text("##gff-version 3\n" + "".join(f"{c}\tbarrnap:0.9\t{t}\t{s}\t{e}\t0\t+\t.\tName={t}\n" for c, t, s, e in rows))
    return str(p)


def read_windows(path):
    return [(l[1:].strip(), None) for l in open(path) if l.startswith(">")], er.read_fasta(path)


def run(tmp_path, h, assembly, gffp, pad=None):
    cmd = [sys.executable, os.path.join(h["BIN"], "extract_rrna_region.py"), "--assembly", assembly, "--gff", gffp, "--out", str(tmp_path / "w.fa")]
    if pad is not None: cmd += ["--pad", str(pad)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return er.read_fasta(str(tmp_path / "w.fa")), r.stdout


def test_read_fasta_takes_first_token_of_header(tmp_path):
    p = tmp_path / "a.fa"; p.write_text(">c1 len=6 cov=3\nACG\nTTT\n>c2\nGG\n")
    assert er.read_fasta(str(p)) == {"c1": "ACGTTT", "c2": "GG"}


def test_overlapping_features_merge_into_one_padded_window(tmp_path, helpers):
    seq = helpers["random_dna"](10000, seed=61)
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"c1": seq})
    # 18S at 4000-5800 and 5.8S at 5900-6050: with pad 100 the windows overlap -> one window 3900-6150
    w, out = run(tmp_path, helpers, asm, gff(tmp_path, [("c1", "rRNA", 4000, 5800), ("c1", "rRNA", 5900, 6050)]), pad=100)
    assert list(w) == ["c1_rRNAregion_0_3900-6150"] and w["c1_rRNAregion_0_3900-6150"] == seq[3900:6150]
    assert "wrote 1 rRNA-operon window" in out


def test_distant_features_give_separate_windows_clipped_to_contig(tmp_path, helpers):
    seq = helpers["random_dna"](5000, seed=62)
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"c1": seq, "c2": seq[:1000]})
    w, _ = run(tmp_path, helpers, asm, gff(tmp_path, [("c1", "rRNA", 50, 200), ("c1", "gene", 4000, 4100), ("c2", "rRNA", 900, 990)]), pad=300)
    assert sorted(w) == ["c1_rRNAregion_0_0-500", "c1_rRNAregion_1_3700-4400", "c2_rRNAregion_0_600-1000"]
    assert w["c2_rRNAregion_0_600-1000"] == seq[600:1000]


def test_non_rrna_rows_and_unknown_contigs_are_ignored(tmp_path, helpers):
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"c1": helpers["random_dna"](2000, seed=63)})
    w, out = run(tmp_path, helpers, asm, gff(tmp_path, [("c1", "CDS", 100, 200), ("ghost", "rRNA", 100, 200)]))
    assert w == {} and "wrote 0" in out


def test_missing_gff_yields_empty_output(tmp_path, helpers):
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"c1": "ACGT" * 100})
    w, _ = run(tmp_path, helpers, asm, str(tmp_path / "absent.gff"))
    assert w == {}
