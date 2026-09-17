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


def test_top_skani_without_ref_name_column_falls_back_to_the_ref_file(tmp_path):
    p = tmp_path / "s.tsv"; p.write_text("Ref_file\tb\tANI\nref1.fa\tq.fa\t99.1\n")
    assert nc.top_skani(str(p)) == (99.1, "ref1.fa")


def run(tmp_path, h, skani_path=None, identify=None, manifest=None):
    cmd = [sys.executable, os.path.join(h["BIN"], "novelty_call.py"), "--sample", "S1", "--out", str(tmp_path / "n.json")]
    if skani_path: cmd += ["--skani", skani_path]
    if manifest: cmd += ["--manifest", manifest]
    if identify is not None:
        p = tmp_path / "id.json"; p.write_text(json.dumps(identify)); cmd += ["--identify-json", str(p)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.load(open(tmp_path / "n.json"))


def test_genome_ani_verdicts(tmp_path, helpers):
    d = run(tmp_path, helpers, skani(tmp_path, [("Aspergillus fumigatus Af293", 99.3)]), {"species": "Aspergillus fumigatus", "its_identity": 99.0})
    assert d["verdict"] == "known_species" and d["basis"] == "genome_ani(skani)" and d["genome_ani"] == 99.3
    assert d["nearest"]["organism"] == "Aspergillus fumigatus Af293" and d["species"] == "Aspergillus fumigatus" and d["flags"] == []
    # a synthetic divergent genome: 88 % to its nearest reference -> nothing close in the set
    d2 = run(tmp_path, helpers, skani(tmp_path, [("X", 88.0)]), {"species": "Aspergillus sp.", "its_identity": 99.0})
    assert d2["verdict"] == "candidate_novel_or_unrepresented" and "no close genome" in d2["note"]      # ANI outranks a good ITS hit
    # 90-95 %: sister lineage of the nearest species unless ITS + a secondary locus vouch for that species
    d3 = run(tmp_path, helpers, skani(tmp_path, [("Aspergillus flavus NRRL3357", 93.5)]), {"species": "Aspergillus flavus", "its_identity": 99.5, "confidence": "medium"})
    assert d3["verdict"] == "candidate_novel_species" and "sister lineage" in d3["note"]
    d4 = run(tmp_path, helpers, skani(tmp_path, [("Aspergillus flavus NRRL3357", 93.5)]), {"species": "Aspergillus flavus", "its_identity": 99.5, "confidence": "high"},
             manifest=None)
    assert d4["verdict"] == "candidate_novel_species"                        # no manifest: the nearest species is unknown, so no rescue
    man = tmp_path / "manifest.tsv"
    man.write_text("accession\torganism\tspecies\tstrain\ttax_id\tcategory\tlevel\tlength\tsource\tfile\nr0\tAspergillus flavus NRRL3357\tAspergillus flavus\tNRRL3357\t5059\treference genome\tComplete Genome\t37000000\tAspergillus\tr0.fa\n")
    d5 = run(tmp_path, helpers, skani(tmp_path, [("Aspergillus flavus NRRL3357", 93.5)]), {"species": "Aspergillus flavus", "its_identity": 99.5, "confidence": "high"}, manifest=str(man))
    assert d5["verdict"] == "known_species" and d5["basis"] == "genome_ani(skani)+ITS+secondary_loci" and d5["nearest"]["species"] == "Aspergillus flavus"
    # identification disagreeing with the nearest genome is flagged, never hidden
    d6 = run(tmp_path, helpers, skani(tmp_path, [("Aspergillus flavus NRRL3357", 98.0)]), {"species": "Aspergillus oryzae", "its_identity": 99.9, "confidence": "high"}, manifest=str(man))
    assert d6["verdict"] == "known_species" and d6["flags"] == ["identification_disagrees:Aspergillus oryzae!=Aspergillus flavus"]


def test_aligned_fraction_floor(tmp_path, helpers):
    p = tmp_path / "s.tsv"
    p.write_text("Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\nr1.fa\tq.fa\t99.0\t3.0\t2.5\tContaminant-like\tq\nr2.fa\tq.fa\t91.0\t60\t55\tSister sp\tq\n")
    rows = nc.read_skani(str(p))
    assert [r["ref_name"] for r in rows] == ["Contaminant-like", "Sister sp"] and rows[0]["af_ok"] is False and rows[1]["af_ok"] is True
    d = run(tmp_path, helpers, str(p), {"species": "X y", "its_identity": 99.0})
    assert d["nearest"]["organism"] == "Sister sp" and d["verdict"] == "candidate_novel_species" and d["n_references_hit"] == 2
    p.write_text("Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\nr1.fa\tq.fa\t99.0\t3.0\t2.5\tContaminant-like\tq\n")
    d2 = run(tmp_path, helpers, str(p), {"species": "X y", "its_identity": 99.0})
    assert d2["verdict"] == "candidate_novel_or_unrepresented" and "aligned fraction" in d2["basis"]


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
