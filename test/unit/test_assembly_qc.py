"""Unit tests for bin/assembly_qc.py — contiguity stats, BUSCO/compleasm parsing, qc_pass."""
import json, os, subprocess, sys
import assembly_qc as aq  # noqa: E402  (sys.path set in conftest)


def test_contig_lengths_and_n50(tmp_path, helpers):
    fa = helpers["write_fasta"](tmp_path / "a.fa", {"c1": "A" * 100, "c2": "C" * 300, "c3": "G" * 50}, width=60)
    lens = aq.contig_lengths(fa)
    assert sorted(lens) == [50, 100, 300]
    assert aq.n50(lens) == 300 and aq.n50([]) == 0 and aq.n50([10, 10, 10, 10]) == 10
    assert aq.contig_lengths(str(tmp_path / "missing.fa")) == [] and aq.contig_lengths(None) == []


def test_busco_complete_prefers_compleasm_then_busco(tmp_path):
    comp = tmp_path / "compleasm"; comp.mkdir()
    (comp / "summary.txt").write_text("## lineage: fungi_odb10\nS:95.50%, 723\nD:2.25%, 17\nF:1.00%, 8\nM:1.25%, 10\nN:758\n")
    bus = tmp_path / "busco"; bus.mkdir()
    (bus / "short_summary.specific.fungi_odb10.busco.txt").write_text("\tC:97.5%[S:96.0%,D:1.5%],F:1.0%,M:1.5%,n:758\n")
    assert aq.busco_complete(str(comp), str(bus)) == 97.75
    assert aq.busco_complete(str(tmp_path / "none"), str(bus)) == 97.5
    assert aq.busco_complete(str(tmp_path / "none"), str(tmp_path / "none2")) is None


def run(tmp_path, h, fasta, **dirs):
    cmd = [sys.executable, os.path.join(h["BIN"], "assembly_qc.py"), "--sample", "S1", "--nuclear", fasta,
           "--out", str(tmp_path / "qc.json"), "--compleasm-dir", dirs.get("compleasm", str(tmp_path / "nocomp")),
           "--busco-dir", dirs.get("busco", str(tmp_path / "nobusco"))]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.load(open(tmp_path / "qc.json")), r.stdout


def test_qc_pass_requires_sequence_and_busco_at_least_80(tmp_path, helpers):
    fa = helpers["write_fasta"](tmp_path / "a.fa", {"c1": "A" * 1000, "c2": "C" * 400})
    bus = tmp_path / "busco"; bus.mkdir()
    (bus / "short_summary.txt").write_text("C:85.0%[S:84.0%,D:1.0%],F:5.0%,M:10.0%,n:758\n")
    d, out = run(tmp_path, helpers, fa, busco=str(bus))
    assert d["assembly_len"] == 1400 and d["n_contigs"] == 2 and d["n50"] == 1000 and d["largest_contig"] == 1000
    assert d["busco_complete"] == 85.0 and d["qc_pass"] is True and d["lineage"] == "fungi_odb10" and "BUSCO 85.0%" in out
    (bus / "short_summary.txt").write_text("C:60.0%[S:59.0%,D:1.0%],F:5.0%,M:35.0%,n:758\n")
    d2, _ = run(tmp_path, helpers, fa, busco=str(bus))
    assert d2["qc_pass"] is False


def test_qc_pass_true_without_completeness_but_false_on_empty_assembly(tmp_path, helpers):
    fa = helpers["write_fasta"](tmp_path / "a.fa", {"c1": "A" * 10})
    d, _ = run(tmp_path, helpers, fa)
    assert d["busco_complete"] is None and d["qc_pass"] is True
    empty = tmp_path / "e.fa"; empty.write_text("")
    d2, _ = run(tmp_path, helpers, str(empty))
    assert d2["assembly_len"] == 0 and d2["qc_pass"] is False and d2["largest_contig"] == 0
