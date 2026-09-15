"""Unit tests for bin/bgc_summary.py — antiSMASH result parsing and the status / note logic."""
import json, os, subprocess, sys
import bgc_summary as bs  # noqa: E402  (sys.path set in conftest)


def as_json(dirpath, regions):
    dirpath.mkdir(exist_ok=True)
    feats = [{"type": "region", "qualifiers": {"product": prods}} for prods in regions]
    feats.append({"type": "CDS", "qualifiers": {"product": ["not a region"]}})
    (dirpath / "sample.json").write_text(json.dumps({"records": [{"features": feats}]}))
    return str(dirpath)


def test_parse_counts_products_and_falls_back_to_region_gbks(tmp_path):
    d = as_json(tmp_path / "as", [["NRPS"], ["T1PKS"], ["NRPS"], ["terpene", "NRPS"]])
    assert bs.parse_antismash(d) == {"NRPS": 3, "T1PKS": 1, "terpene": 1}
    e = tmp_path / "as2"; e.mkdir()
    for i in range(3):
        (e / f"c1.region00{i}.gbk").write_text("LOCUS x\n")
    assert bs.parse_antismash(str(e)) == {"region": 3}
    assert bs.parse_antismash(str(tmp_path / "empty")) == {}
    (tmp_path / "as3").mkdir(); (tmp_path / "as3" / "bad.json").write_text("{nope")
    assert bs.parse_antismash(str(tmp_path / "as3")) == {}


def run(tmp_path, h, as_dir, status=None, log=None):
    cmd = [sys.executable, os.path.join(h["BIN"], "bgc_summary.py"), "--sample", "S1", "--as-dir", as_dir, "--out", str(tmp_path / "bgc.json")]
    if status: cmd += ["--status", status]
    if log: cmd += ["--log", log]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.load(open(tmp_path / "bgc.json")), r.stdout


def test_ok_run_reports_counts(tmp_path, helpers):
    d, out = run(tmp_path, helpers, as_json(tmp_path / "as", [["NRPS"], ["T1PKS"], ["NRPS"]]), status="ok")
    assert d["status"] == "ok" and d["n_clusters"] == 3 and d["by_type"] == {"NRPS": 2, "T1PKS": 1}
    assert d["clusters"] == [{"type": "NRPS", "count": 2}, {"type": "T1PKS", "count": 1}] and "note" not in d
    assert "3 BGCs" in out


def test_ok_status_without_result_json_becomes_failed(tmp_path, helpers):
    e = tmp_path / "as"; e.mkdir()
    d, _ = run(tmp_path, helpers, str(e), status="ok")
    assert d["status"] == "failed" and d["n_clusters"] is None and d["note"] == "antiSMASH failed"


def test_failed_status_records_log_tail(tmp_path, helpers):
    log = tmp_path / "as.log"; log.write_text("\n".join(f"line {i}" for i in range(20)) + "\n")
    d, _ = run(tmp_path, helpers, str(tmp_path / "none"), status="failed", log=str(log))
    assert d["status"] == "failed" and d["n_clusters"] is None
    assert d["note"].startswith("antiSMASH failed: line 12") and d["note"].endswith("line 19")


def test_skipped_status_and_inferred_status(tmp_path, helpers):
    d, _ = run(tmp_path, helpers, str(tmp_path / "none"), status="skipped")
    assert d["status"] == "skipped" and d["n_clusters"] is None
    d2, _ = run(tmp_path, helpers, as_json(tmp_path / "as", [["NRPS"]]))     # no --status: inferred from the JSON presence
    assert d2["status"] == "ok" and d2["n_clusters"] == 1
