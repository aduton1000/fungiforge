"""Unit tests for bin/validate_run.py — the expected-results comparison (W1.1)."""
import json, os, subprocess, sys
import validate_run as vr  # noqa: E402  (sys.path set in conftest)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def make_results(tmp_path, sample="S1", **over):
    """A minimal results tree: stage JSONs, master row, provenance."""
    stages = {
        "assembly_qc": {"stage": "assembly_qc", "assembly_len": 29900000, "n_contigs": 43, "n50": 1900000, "busco_complete": 99.7, "qc_pass": True},
        "identify": {"stage": "identify", "species": "Aspergillus fumigatus", "confidence": "high", "its_present": True, "its_identity": 100.0},
        "resistance": {"stage": "resistance", "summary": {"n_known_mutations": 0, "resistant_drug_classes": []}, "cyp51A_TR": {"tr_type": "none", "cyp51A_located": True}},
        "bgc": {"stage": "bgc", "status": "ok", "n_clusters": 50},
        "decontam": {"stage": "decontam", "verdict": "fungal", "dropped_bp_pct": 0.1},
        "annotate": {"stage": "annotate", "n_proteins": 9630},
        "polish": {"stage": "polish", "mode": "hybrid"},
        "novelty": {"stage": "novelty", "verdict": "known_species"},
    }
    for k, v in over.items():
        st, key = k.split("__")
        if st != "master":
            stages[st][key] = v
    dirs = {"assembly_qc": "05_assembly_qc", "identify": "08_identify", "resistance": "09_resistance", "bgc": "11_bgc",
            "decontam": "04_decontam", "annotate": "07_annotate", "polish": "03_polish", "novelty": "12_novelty"}
    for st, doc in stages.items():
        d = tmp_path / "results" / sample / dirs[st]; d.mkdir(parents=True, exist_ok=True)
        (d / f"{sample}.{st}.json").write_text(json.dumps(doc))
    ms = tmp_path / "results" / "04_summary"; ms.mkdir(exist_ok=True)
    (ms / f"{sample}.master.tsv").write_text("sample\tresistant_classes\tstages_failed\n" + f"{sample}\tnone\t{over.get('master__stages_failed', 'none')}\n")
    pi = tmp_path / "results" / "pipeline_info"; pi.mkdir(exist_ok=True)
    (pi / "provenance.json").write_text(json.dumps({"n_samples": 1}))
    return str(tmp_path / "results")


def cea10_expected():
    with open(os.path.join(ROOT, "test", "expected", "AfumCEA10.json")) as fh:
        e = json.load(fh)
    e["sample"] = "S1"
    return e


def test_evaluate_rules():
    assert vr.evaluate({"equals": "none"}, "none", True) == ("PASS", "")
    assert vr.evaluate({"equals": "none"}, "bgc:failed", True)[0] == "FAIL"
    assert vr.evaluate({"equals": True}, "True", True)[0] == "PASS"          # master TSV strings compare case-insensitively
    assert vr.evaluate({"one_of": ["fungal", "likely_fungal"]}, "mixed", True)[0] == "FAIL"
    assert vr.evaluate({"regex": "^Aspergillus "}, "Aspergillus flavus", True)[0] == "PASS"
    assert vr.evaluate({"min": 98}, 99.7, True)[0] == "PASS" and vr.evaluate({"min": 98}, "97.9", True)[0] == "FAIL"
    assert vr.evaluate({"max": 60}, 43, True)[0] == "PASS" and vr.evaluate({"max": 60}, 61, True)[0] == "FAIL"
    assert vr.evaluate({"expected": 100, "rel_tol": 0.05}, 104, True)[0] == "PASS"
    assert vr.evaluate({"expected": 100, "rel_tol": 0.05}, 106, True)[0] == "FAIL"
    assert vr.evaluate({"expected": 100, "abs_tol": 10}, 109, True)[0] == "PASS"
    assert vr.evaluate({"min": 1}, "not a number", True)[0] == "FAIL"
    assert vr.evaluate({"equals": 1}, None, False) == ("MISSING", "value not found")
    assert vr.evaluate({"equals": 1, "optional": True}, None, False)[0] == "PASS"


def test_reference_run_passes_every_check(tmp_path):
    res = make_results(tmp_path)
    rep = vr.validate(res, cea10_expected())
    assert rep["passed"], [c for c in rep["checks"] if c["status"] != "PASS"]
    assert rep["n_checks"] == len(cea10_expected()["checks"])


def test_deviations_are_reported_by_name(tmp_path):
    res = make_results(tmp_path, assembly_qc__assembly_len=31000000, bgc__status="failed", identify__its_present=False, master__stages_failed="bgc:failed")
    rep = vr.validate(res, cea10_expected())
    failed = {c["check"]: c for c in rep["checks"] if c["status"] != "PASS"}
    assert set(failed) == {"assembly_len", "bgc_status", "its_present", "stages_all_ok"}
    assert "±" in failed["assembly_len"]["detail"]
    assert rep["n_fail"] == 4 and not rep["passed"]


def test_missing_stage_json_is_missing_not_crash(tmp_path):
    res = make_results(tmp_path)
    os.remove(os.path.join(res, "S1", "11_bgc", "S1.bgc.json"))
    rep = vr.validate(res, cea10_expected())
    st = {c["check"]: c["status"] for c in rep["checks"]}
    assert st["bgc_status"] == "MISSING" and st["n_bgc"] == "MISSING" and st["provenance_samples"] == "PASS"


def test_cli_writes_table_and_exit_status(tmp_path):
    res = make_results(tmp_path)
    exp = tmp_path / "exp.json"; exp.write_text(json.dumps(cea10_expected()))
    r = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "validate_run.py"), "--results", res, "--expected", str(exp),
                        "--out", str(tmp_path / "v.tsv"), "--json", str(tmp_path / "v.json")], capture_output=True, text=True)
    assert r.returncode == 0 and "checks passed" in r.stdout
    rows = [l.split("\t") for l in open(tmp_path / "v.tsv").read().splitlines()]
    assert rows[0][:3] == ["sample", "check", "from"] and all(row[5] == "PASS" for row in rows[1:])
    assert json.load(open(tmp_path / "v.json"))["passed"] is True
    res2 = make_results(tmp_path / "b", identify__species="Aspergillus flavus")
    r2 = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "validate_run.py"), "--results", res2, "--expected", str(exp)], capture_output=True, text=True)
    assert r2.returncode == 1 and "FAIL    species" in r2.stdout


def test_shipped_expected_files_are_well_formed():
    for name in ("AfumCEA10.json", "ASSARM-PHI-DF-005.json"):
        with open(os.path.join(ROOT, "test", "expected", name)) as fh:
            e = json.load(fh)
        assert e["sample"] and e["checks"]
        for k, rule in e["checks"].items():
            assert "from" in rule and "." in rule["from"], k
            assert any(r in rule for r in ("equals", "one_of", "regex", "min", "max", "expected")), k
