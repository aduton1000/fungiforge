"""Unit tests for bin/id_classify.py — UNITE best-hit parsing, ITS thresholds, sourmash concordance."""
import json, os, subprocess, sys
import id_classify as idc  # noqa: E402  (sys.path set in conftest)

UNITE = "Aspergillus_fumigatus|KC411932|SH1136727.09FU|refs|k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Aspergillus;s__Aspergillus_fumigatus"
UNITE_SP = "Aspergillus_sp|XX000001|SH0000001.09FU|reps|k__Fungi;p__Ascomycota;g__Aspergillus;s__Aspergillus_sp"


def b6(tmp_path, rows):
    p = tmp_path / "unite.b6"
    p.write_text("".join(f"its1\t{t}\t{pid}\t500\t2\t0\t1\t500\t1\t500\t0\t900\n" for t, pid in rows))
    return str(p)


def test_parse_unite_b6_takes_best_identity_and_splits_taxonomy(tmp_path):
    hit = idc.parse_unite_b6(b6(tmp_path, [(UNITE_SP, 97.0), (UNITE, 99.6), (UNITE_SP, 95.0)]))
    assert hit == {"pident": 99.6, "species": "Aspergillus fumigatus", "genus": "Aspergillus", "sh": "SH1136727.09FU"}


def test_parse_unite_b6_drops_sp_and_unidentified_species(tmp_path):
    hit = idc.parse_unite_b6(b6(tmp_path, [(UNITE_SP, 99.0)]))
    assert hit["species"] is None and hit["genus"] == "Aspergillus"
    assert idc.parse_unite_b6(None) is None and idc.parse_unite_b6(str(tmp_path / "none.b6")) is None


def test_top_gather_extracts_binomial(tmp_path):
    g = tmp_path / "g.csv"
    g.write_text("intersect_bp,f_orig_query,name\n1000,0.9,\"GCF_000002655.1 Aspergillus fumigatus Af293 strain=Af293\"\n")
    assert idc.top_gather(str(g))["species"] == "Aspergillus fumigatus"
    g.write_text("intersect_bp,name\n")
    assert idc.top_gather(str(g)) is None


def run(tmp_path, h, its_text=">its1\nACGT\n", unite=None, gather=None):
    its = tmp_path / "its.fa"; its.write_text(its_text)
    cmd = [sys.executable, os.path.join(h["BIN"], "id_classify.py"), "--sample", "S1", "--its", str(its),
           "--out-species", str(tmp_path / "sp.txt"), "--out-json", str(tmp_path / "id.json")]
    if unite: cmd += ["--unite-b6", unite]
    if gather: cmd += ["--gather", gather]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.load(open(tmp_path / "id.json")), open(tmp_path / "sp.txt").read().strip().split("\t")


def test_species_level_call_at_or_above_98_5(tmp_path, helpers):
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 98.5)]))
    assert sp == ["Aspergillus fumigatus", "high"] and doc["method"] == "ITS/UNITE" and doc["its_identity"] == 98.5
    assert doc["its_present"] is True and doc["stage"] == "identify"


def test_medium_confidence_between_94_and_98_5(tmp_path, helpers):
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 96.0)]))
    assert sp == ["Aspergillus fumigatus", "medium"] and "below-species-threshold" in doc["method"]


def test_genus_only_below_94_or_sp_hit(tmp_path, helpers):
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 90.0)]))
    assert sp == ["Aspergillus sp.", "low"] and "genus-only" in doc["method"]
    doc2, sp2 = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE_SP, 99.0)]))
    assert sp2 == ["Aspergillus sp.", "low"]


def test_unknown_without_hits_and_empty_its(tmp_path, helpers):
    doc, sp = run(tmp_path, helpers, its_text="")
    assert sp == ["unknown", "none"] and doc["method"] == "no_reference" and doc["its_present"] is False and doc["its_identity"] is None


def test_gather_concordance_raises_confidence_and_fallback(tmp_path, helpers):
    g = tmp_path / "g.csv"
    g.write_text("intersect_bp,name\n1000,\"GCF_1 Aspergillus fumigatus Af293\"\n")
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 96.0)]), gather=str(g))
    assert sp == ["Aspergillus fumigatus", "high"] and doc["method"].endswith("+genome_ani(concordant/GCPSR)")
    doc2, sp2 = run(tmp_path, helpers, its_text="", gather=str(g))
    assert sp2 == ["Aspergillus fumigatus", "medium"] and doc2["method"] == "genome_ani(sourmash)"
    g.write_text("intersect_bp,name\n1000,\"GCF_2 Candida albicans SC5314\"\n")
    doc3, sp3 = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 99.0)]), gather=str(g))
    assert sp3 == ["Aspergillus fumigatus", "high"] and "concordant" not in doc3["method"]   # discordant genome hit is not merged
