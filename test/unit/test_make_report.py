"""Unit tests for bin/make_report.py — master-table row construction and the status columns."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "bin")))
import make_report as m  # noqa: E402


def row(S):
    return m.build_row("X", "AIR", "F1", "Dry", S)


def test_columns_are_append_only_and_end_with_status_columns():
    assert m.MASTER_COLS[:4] == ["sample", "compartment", "facility", "season"]
    assert m.MASTER_COLS[-5:] == ["sample_verdict", "contam_removed_pct", "top_taxon", "stages_failed", "gate"]
    assert len(m.MASTER_COLS) == len(set(m.MASTER_COLS))


def test_row_has_every_column():
    r = row({})
    assert set(r) == set(m.MASTER_COLS)


def test_stages_failed_none_when_all_ok():
    S = {"annotate": {"status": "ok"}, "bgc": {"status": "ok", "n_clusters": 3}}
    assert row(S)["stages_failed"] == "none"


def test_stages_failed_lists_non_ok_stages_sorted():
    S = {"bgc": {"status": "failed"}, "polish": {"status": "partial"}, "novelty": {"status": "skipped"}, "annotate": {"status": "ok"}}
    assert row(S)["stages_failed"] == "bgc:failed;novelty:skipped;polish:partial"


def test_stages_without_status_are_not_flagged():
    S = {"readqc": {"illumina": True}}
    assert row(S)["stages_failed"] == "none"


def test_bgc_failed_gives_na_not_zero():
    assert row({"bgc": {"status": "failed", "n_clusters": None}})["n_bgc"] == "NA"
    assert row({"bgc": {"status": "ok", "n_clusters": 51}})["n_bgc"] == 51


def test_decontam_columns():
    S = {"decontam": {"verdict": "fungal", "dropped_bp_pct": 0.04,
                      "top_species": [{"name": "Aspergillus flavus", "pct": 51.3}]}}
    r = row(S)
    assert r["sample_verdict"] == "fungal" and r["contam_removed_pct"] == 0.04
    assert r["top_taxon"] == "Aspergillus flavus (51.3%)"


# ---------------------------------------------------------------- CLI end to end (W0.3)
import json, subprocess, sys  # noqa: E402


def test_cli_merges_stage_jsons_into_html_and_master_row(tmp_path):
    stages = {
        "identify": {"stage": "identify", "species": "Aspergillus fumigatus", "confidence": "high", "method": "ITS/UNITE", "status": "ok"},
        "assembly_qc": {"stage": "assembly_qc", "qc_pass": True, "busco_complete": 99.1, "lineage": "fungi_odb10", "assembly_len": 29000000, "n_contigs": 43, "status": "ok"},
        "resistance": {"stage": "resistance", "status": "ok", "polish_mode": "hybrid",
                       "calls": [{"gene": "cyp51A", "drug_class": "azole", "change": "L98H", "known": True, "confidence": "high"}],
                       "cyp51A_TR": {"tr_type": "TR34"}, "summary": {"resistant_drug_classes": ["azole"], "n_known_mutations": 1}},
        "bgc": {"stage": "bgc", "status": "failed", "n_clusters": None},
        "decontam": {"stage": "decontam", "status": "ok", "verdict": "fungal", "dropped_bp_pct": 0.5, "top_species": [{"name": "Aspergillus fumigatus", "pct": 30.0}]},
    }
    for k, v in stages.items():
        (tmp_path / f"S1.{k}.json").write_text(json.dumps(v))
    (tmp_path / "S1.broken.json").write_text("{not json")
    r = subprocess.run([sys.executable, os.path.join(HERE, "..", "..", "bin", "make_report.py"), "--sample", "S1", "--compartment", "AIR",
                        "--facility", "F1", "--season", "Dry", "--jsons", str(tmp_path / "S1.*.json"),
                        "--html", str(tmp_path / "S1.html"), "--master", str(tmp_path / "S1.master.tsv")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    header, row = [l.rstrip("\n").split("\t") for l in open(tmp_path / "S1.master.tsv")]
    assert header == m.MASTER_COLS
    d = dict(zip(header, row))
    assert d["sample"] == "S1" and d["species"] == "Aspergillus fumigatus" and d["qc_pass"] == "True"
    assert d["resistant_classes"] == "azole" and d["cyp51A_TR"] == "TR34" and d["n_bgc"] == "NA" and d["stages_failed"] == "bgc:failed"
    assert d["top_taxon"] == "Aspergillus fumigatus (30.0%)" and d["sample_verdict"] == "fungal" and d["contam_removed_pct"] == "0.5"
    html_text = open(tmp_path / "S1.html").read()
    assert "<title>FungiForge · S1</title>" in html_text and "L98H" in html_text and "Aspergillus fumigatus" in html_text
    assert "5 stages" in r.stdout          # the unreadable JSON was skipped, not fatal



def test_gate_column_reports_skipped_reason_and_stages_failed():
    S = {"decontam": {"status": "ok", "verdict": "non_fungal"},
         "gate": {"stage": "gate", "status": "skipped", "reason": "verdict:non_fungal"}}
    r = row(S)
    assert r["gate"] == "skipped(verdict:non_fungal)" and r["stages_failed"] == "gate:skipped"
    assert r["species"] == "NA" and r["n_bgc"] == "NA"
    assert row({"decontam": {"status": "ok", "verdict": "fungal"}})["gate"] == "pass"
