"""Unit tests for bin/busco_lineage.py — species-aware BUSCO summary (stage 08b)."""
import json, os, subprocess, sys
import busco_lineage as bl  # noqa: E402  (sys.path set in conftest)


def test_breakdown_from_compleasm_then_busco(tmp_path):
    comp = tmp_path / "compleasm"; comp.mkdir()
    (comp / "summary.txt").write_text("## lineage: eurotiales_odb10\nS:97.50%, 4100\nD:1.25%, 52\nF:0.50%, 21\nM:0.75%, 31\nN:4204\n")
    assert bl.busco_breakdown(str(comp), None) == {"single": 97.5, "duplicated": 1.25, "fragmented": 0.5, "missing": 0.75}
    bus = tmp_path / "busco"; bus.mkdir()
    (bus / "short_summary.specific.eurotiales_odb10.busco.txt").write_text("\tC:98.9%[S:98.0%,D:0.9%],F:0.4%,M:0.7%,n:4191\n")
    assert bl.busco_breakdown(str(tmp_path / "none"), str(bus)) == {"single": 98.0, "duplicated": 0.9, "fragmented": 0.4, "missing": 0.7}
    assert bl.busco_breakdown(None, None) == {}


def run(tmp_path, h, *args):
    out = tmp_path / "bl.json"
    r = subprocess.run([sys.executable, os.path.join(h["BIN"], "busco_lineage.py"), "--sample", "S1", "--out", str(out), *args],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.load(open(out)), r.stdout


def test_cli_skipped_and_scored(tmp_path, helpers):
    doc, out = run(tmp_path, helpers, "--lineage", "fungi_odb10", "--qc-lineage", "fungi_odb10", "--skipped", "same_as_qc")
    assert doc["stage"] == "busco_lineage" and doc["busco_complete"] is None and doc["skipped"] == "same_as_qc" and "skipped" in out
    comp = tmp_path / "compleasm"; comp.mkdir()
    (comp / "summary.txt").write_text("S:97.50%, 4100\nD:1.25%, 52\nF:0.50%, 21\nM:0.75%, 31\n")
    doc2, _ = run(tmp_path, helpers, "--lineage", "eurotiales_odb10", "--compleasm-dir", str(comp), "--busco-dir", str(tmp_path / "nb"))
    assert doc2["lineage"] == "eurotiales_odb10" and doc2["busco_complete"] == 98.75 and doc2["single"] == 97.5 and doc2["skipped"] is None
    doc3, _ = run(tmp_path, helpers, "--lineage", "eurotiales_odb10", "--compleasm-dir", str(tmp_path / "x"), "--busco-dir", str(tmp_path / "y"))
    assert doc3["busco_complete"] is None and "no summary" in doc3["note"]
