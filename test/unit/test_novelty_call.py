"""Unit tests for bin/novelty_call.py — skani parsing and the ANI / ITS / undetermined verdicts."""
import json, os, subprocess, sys
import novelty_call as nc  # noqa: E402  (sys.path set in conftest)


def skani(tmp_path, rows, header="Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name"):
    p = tmp_path / "skani.tsv"
    p.write_text(header + "\n" + "".join(f"r{i}.fa\tq.fa\t{ani}\t90\t90\t{name}\tq\n" for i, (name, ani) in enumerate(rows)))
    return str(p)


def test_top_skani_picks_highest_ani(tmp_path):
    assert nc.top_skani(skani(tmp_path, [("A", 91.2), ("B", 97.8), ("C", 80.0)])) == (97.8, "B")
    assert nc.top_skani(str(tmp_path / "absent.tsv")) is None and nc.top_skani(None) is None
    assert nc.top_skani(skani(tmp_path, [])) is None


def test_top_skani_without_ref_name_column_falls_back_to_second_column(tmp_path):
    p = tmp_path / "s.tsv"; p.write_text("a\tb\tANI\nref1.fa\tq.fa\t99.1\n")
    assert nc.top_skani(str(p)) == (99.1, "q.fa")


def run(tmp_path, h, skani_path=None, identify=None):
    cmd = [sys.executable, os.path.join(h["BIN"], "novelty_call.py"), "--sample", "S1", "--out", str(tmp_path / "n.json")]
    if skani_path: cmd += ["--skani", skani_path]
    if identify is not None:
        p = tmp_path / "id.json"; p.write_text(json.dumps(identify)); cmd += ["--identify-json", str(p)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.load(open(tmp_path / "n.json"))


def test_genome_ani_verdicts(tmp_path, helpers):
    d = run(tmp_path, helpers, skani(tmp_path, [("Aspergillus fumigatus Af293", 99.3)]), {"species": "Aspergillus fumigatus", "its_identity": 99.0})
    assert d["verdict"] == "known_species" and d["basis"] == "genome_ani(skani)" and d["genome_ani"] == 99.3
    assert d["closest_reference"] == "Aspergillus fumigatus Af293" and d["species"] == "Aspergillus fumigatus"
    d2 = run(tmp_path, helpers, skani(tmp_path, [("X", 88.0)]), {"species": "Aspergillus sp.", "its_identity": 99.0})
    assert d2["verdict"] == "candidate_novel_species" and "ANI < 95%" in d2["note"]      # ANI outranks a good ITS hit


def test_its_fallback_verdicts(tmp_path, helpers):
    d = run(tmp_path, helpers, identify={"species": "Aspergillus fumigatus", "its_identity": 99.2})
    assert d["verdict"] == "known_species" and d["basis"] == "ITS_distance(UNITE)" and d["genome_ani"] is None
    d2 = run(tmp_path, helpers, identify={"species": "Aspergillus sp.", "its_identity": 96.0})
    assert d2["verdict"] == "candidate_novel_species" and "< 98.5%" in d2["note"]


def test_undetermined_without_any_signal(tmp_path, helpers):
    d = run(tmp_path, helpers)
    assert d["verdict"] == "undetermined" and d["basis"] == "none" and d["stage"] == "novelty"
    d2 = run(tmp_path, helpers, identify={"species": "unknown", "its_identity": None})
    assert d2["verdict"] == "undetermined"
