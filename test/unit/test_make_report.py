"""Unit tests for bin/make_report.py — master-table row construction and the status columns."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "bin")))
import make_report as m  # noqa: E402


def row(S):
    return m.build_row("X", "AIR", "F1", "Dry", S)


def test_columns_are_append_only_and_end_with_status_columns():
    assert m.MASTER_COLS[:4] == ["sample", "compartment", "facility", "season"]
    assert m.MASTER_COLS[-4:] == ["sample_verdict", "contam_removed_pct", "top_taxon", "stages_failed"]
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
