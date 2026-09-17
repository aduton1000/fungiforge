"""Unit tests for the small stage helpers: bin/te_summary.py, bin/mobile_merge.py."""
import json, os, subprocess, sys


def run(h, script, *args):
    r = subprocess.run([sys.executable, os.path.join(h["BIN"], script), *args], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


# ---------------------------------------------------------------- te_summary
def test_te_summary_counts_families_by_class(tmp_path, helpers):
    telib = helpers["write_fasta"](tmp_path / "telib.fa", {"rnd-1_family-1#LTR/Gypsy": "A", "rnd-1_family-2#LTR/Gypsy": "C",
                                                            "rnd-2_family-3#DNA/TcMar": "G", "rnd-3_family-4": "T"})
    out = run(helpers, "te_summary.py", "--telib", telib, "--out", str(tmp_path / "te.json"))
    d = json.load(open(tmp_path / "te.json"))
    assert d["n_families"] == 4 and d["by_class"] == {"LTR/Gypsy": 2, "DNA/TcMar": 1, "Unknown": 1} and "4 TE families" in out


def test_te_summary_missing_library_gives_zero(tmp_path, helpers):
    run(helpers, "te_summary.py", "--telib", str(tmp_path / "absent.fa"), "--out", str(tmp_path / "te.json"))
    assert json.load(open(tmp_path / "te.json"))["n_families"] == 0


# ---------------------------------------------------------------- mobile_merge
def test_mobile_merge_combines_te_and_genomad(tmp_path, helpers):
    te = tmp_path / "te.json"; te.write_text(json.dumps({"n_families": 7, "by_class": {"LTR/Copia": 7}}))
    g = tmp_path / "genomad_out" / "x_summary"; g.mkdir(parents=True)
    (g / "x_virus_summary.tsv").write_text("seq_name\tlength\n v1\t100\n v2\t200\n")
    out = run(helpers, "mobile_merge.py", "--sample", "S1", "--te", str(te), "--genomad", str(tmp_path / "genomad_out"), "--out", str(tmp_path / "m.json"))
    d = json.load(open(tmp_path / "m.json"))
    assert d["n_te_families"] == 7 and d["te_by_class"] == {"LTR/Copia": 7} and d["n_mycovirus"] == 2
    assert d["stage"] == "mobile" and "DNA-only" in d["note"] and "2 viral contigs" in out


def test_mobile_merge_degrades_without_inputs(tmp_path, helpers):
    run(helpers, "mobile_merge.py", "--sample", "S1", "--out", str(tmp_path / "m.json"))
    d = json.load(open(tmp_path / "m.json"))
    assert d["n_te_families"] == "NA" and d["te_by_class"] == {} and d["n_mycovirus"] == 0

# extras (stage 13) is tested in test_extras.py since W2.6 (tool-backed, no keyword tallies)
