"""Unit tests for the stage status contract (bin/ff_status.py + bin/ff_status.sh)."""
import json, os, subprocess, sys, textwrap
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.abspath(os.path.join(HERE, "..", "..", "bin"))
sys.path.insert(0, BIN)
import ff_status  # noqa: E402


# ---------------------------------------------------------------- compute_status
def test_status_ok_when_all_tools_exit_zero():
    tools = {"a": {"exit": 0, "optional": False, "seconds": 1}, "b": {"exit": 0, "optional": True, "seconds": 1}}
    assert ff_status.compute_status(tools) == ("ok", "")


def test_status_partial_when_only_optional_tool_fails():
    tools = {"a": {"exit": 0, "optional": False, "seconds": 1}, "b": {"exit": 2, "optional": True, "seconds": 1}}
    status, note = ff_status.compute_status(tools)
    assert status == "partial" and "b (exit 2)" in note


def test_status_failed_when_required_tool_fails_even_if_optional_ok():
    tools = {"a": {"exit": 1, "optional": False, "seconds": 1}, "b": {"exit": 0, "optional": True, "seconds": 1}}
    status, note = ff_status.compute_status(tools)
    assert status == "failed" and "a (exit 1)" in note


def test_status_skipped_overrides_everything():
    tools = {"a": {"exit": 1, "optional": False, "seconds": 1}}
    assert ff_status.compute_status(tools, stage_skipped="disabled")[0] == "skipped"


def test_no_tools_is_ok():
    assert ff_status.compute_status({}) == ("ok", "")


# ---------------------------------------------------------------- finalize CLI
def run_finalize(tmp_path, tools_lines, best_effort=False, json_doc=None, skips=None):
    tools = tmp_path / "t.tsv"; tools.write_text("".join(tools_lines))
    sk = tmp_path / "s.tsv"; sk.write_text("".join(skips or []))
    js = tmp_path / "stage.json"
    if json_doc is not None:
        js.write_text(json.dumps(json_doc))
    cmd = [sys.executable, os.path.join(BIN, "ff_status.py"), "finalize", "--stage", "demo", "--sample", "S1",
           "--json", str(js), "--tools", str(tools), "--skips", str(sk)] + (["--best-effort"] if best_effort else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode, json.load(open(js)), r.stdout


def test_finalize_merges_into_existing_json_and_exits_zero_when_ok(tmp_path):
    rc, doc, out = run_finalize(tmp_path, ["flye\t0\t0\t12\n"], json_doc={"sample": "S1", "stage": "demo", "n": 3})
    assert rc == 0 and doc["status"] == "ok" and doc["n"] == 3
    assert doc["tools"]["flye"] == {"exit": 0, "optional": False, "seconds": 12}


def test_finalize_creates_json_when_stage_writer_never_ran(tmp_path):
    rc, doc, _ = run_finalize(tmp_path, ["flye\t137\t0\t5\n"])
    assert rc == 1 and doc["status"] == "failed" and doc["sample"] == "S1" and doc["stage"] == "demo"
    assert "flye (exit 137)" in doc["note"]


def test_finalize_best_effort_records_failure_but_exits_zero(tmp_path):
    rc, doc, _ = run_finalize(tmp_path, ["antismash\t1\t0\t5\n"], best_effort=True)
    assert rc == 0 and doc["status"] == "failed"


def test_finalize_partial_exits_zero_and_lists_fallback(tmp_path):
    rc, doc, _ = run_finalize(tmp_path, ["bwa\t0\t1\t1\n", "polypolish\t3\t1\t2\n"])
    assert rc == 0 and doc["status"] == "partial" and "polypolish (exit 3)" in doc["note"]


def test_finalize_records_skips_without_changing_status(tmp_path):
    rc, doc, _ = run_finalize(tmp_path, ["kraken2\t0\t0\t4\n"], skips=["genomad\tdatabase not staged\n"])
    assert rc == 0 and doc["status"] == "ok" and doc["skipped_tools"] == {"genomad": "database not staged"}


def test_finalize_unreadable_json_becomes_failed(tmp_path):
    js = tmp_path / "stage.json"; js.write_text("{not json")
    tools = tmp_path / "t.tsv"; tools.write_text("x\t0\t0\t1\n")
    r = subprocess.run([sys.executable, os.path.join(BIN, "ff_status.py"), "finalize", "--stage", "demo", "--sample", "S1",
                        "--json", str(js), "--tools", str(tools)], capture_output=True, text=True)
    doc = json.load(open(js))
    assert r.returncode == 1 and doc["status"] == "failed" and "unreadable" in doc["note"]


# ---------------------------------------------------------------- bash helper end to end
def run_bash(tmp_path, body, best_effort=""):
    script = tmp_path / "t.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash -ue
        source "{BIN}/ff_status.sh"
        ff_init demo S1 stage.json {best_effort}
        {body}
        ff_finalize
        """))
    r = subprocess.run(["bash", str(script)], cwd=tmp_path, capture_output=True, text=True)
    doc = json.load(open(tmp_path / "stage.json")) if (tmp_path / "stage.json").exists() else None
    return r.returncode, doc, r.stderr


def test_bash_ok_path_records_tool_and_exit_zero(tmp_path):
    rc, doc, _ = run_bash(tmp_path, "ff_run echo_ok -- echo hello\necho '{\"stage\":\"demo\",\"sample\":\"S1\",\"x\":1}' > stage.json")
    assert rc == 0 and doc["status"] == "ok" and doc["tools"]["echo_ok"]["exit"] == 0 and doc["x"] == 1


def test_bash_required_failure_exits_nonzero_and_stops_script(tmp_path):
    rc, doc, err = run_bash(tmp_path, "ff_run boom -- false\necho SHOULD_NOT_RUN > marker.txt")
    assert rc == 1 and doc["status"] == "failed" and doc["tools"]["boom"]["exit"] == 1
    assert not (tmp_path / "marker.txt").exists() and "required step" in err


def test_bash_required_failure_in_best_effort_stage_exits_zero(tmp_path):
    rc, doc, _ = run_bash(tmp_path, "ff_run boom -- false", best_effort="--best-effort")
    assert rc == 0 and doc["status"] == "failed"


def test_bash_optional_failure_continues_and_marks_partial(tmp_path):
    rc, doc, _ = run_bash(tmp_path, "ff_run opt --optional -- false\nif [ \"$FF_RC\" -ne 0 ]; then echo B > m; else echo A > m; fi\nff_run after -- true")
    assert rc == 0 and doc["status"] == "partial" and (tmp_path / "m").read_text().strip() == "B"
    assert doc["tools"]["after"]["exit"] == 0


def test_bash_bare_optional_failure_does_not_trip_errexit(tmp_path):
    # Nextflow runs task scripts with `bash -ue`: a bare `ff_run --optional` statement whose
    # tool fails must NOT abort the script (this was a real bug caught in-container).
    rc, doc, _ = run_bash(tmp_path, "ff_run opt --optional -- false\necho reached > m")
    assert rc == 0 and doc["status"] == "partial" and (tmp_path / "m").exists()


def test_bash_pipeline_exit_code_is_captured_with_pipefail(tmp_path):
    rc, doc, _ = run_bash(tmp_path, "ff_run pipe --optional -- bash -o pipefail -c 'false | cat'")
    assert rc == 0 and doc["tools"]["pipe"]["exit"] != 0 and doc["status"] == "partial"


def test_bash_skip_is_recorded(tmp_path):
    """A stage whose only outcome is a skip reports `skipped`: nothing ran, so nothing succeeded."""
    rc, doc, _ = run_bash(tmp_path, "ff_skip genomad 'database not staged'")
    assert rc == 0 and doc["status"] == "skipped"
    assert doc["skipped_tools"]["genomad"] == "database not staged"


def test_bash_skip_alongside_real_work_stays_ok(tmp_path):
    """The stage did its job; the skipped extra is recorded without downgrading the status."""
    rc, doc, _ = run_bash(tmp_path, "ff_skip genemark 'no licence key'\nff_run echo -- echo hi")
    assert rc == 0 and doc["status"] == "ok"
    assert doc["skipped_tools"]["genemark"] == "no licence key"


# ---------------------------------------------------------------- ff_version (W0.2)
def test_bash_version_is_recorded_and_copied_onto_matching_tool(tmp_path):
    body = "ff_version echo -- echo 'demo tool 1.2.3 (build 7)'\nff_run echo -- echo hi"
    rc, doc, _ = run_bash(tmp_path, body)
    assert rc == 0 and doc["versions"] == {"echo": "demo tool 1.2.3 (build 7)"}
    assert doc["tools"]["echo"]["version"] == "demo tool 1.2.3 (build 7)"


def test_bash_version_picks_first_version_looking_line(tmp_path):
    rc, doc, _ = run_bash(tmp_path, "ff_version t -- printf 'Program: tool\\nUsage: tool [opts]\\nVersion: 0.7.19-r1273\\n'")
    assert doc["versions"]["t"] == "Version: 0.7.19-r1273"


def test_bash_version_falls_back_to_first_line_then_unknown(tmp_path):
    rc, doc, _ = run_bash(tmp_path, "ff_version a -- printf 'no digits here\\n'\nff_version b -- true")
    assert doc["versions"] == {"a": "no digits here", "b": "unknown"}


def test_bash_version_missing_command_and_failing_command_never_abort(tmp_path):
    body = "ff_version gone -- definitely_not_a_command_xyz --version\nff_version bad -- bash -c 'echo v9.9 >&2; exit 3'\necho reached > m"
    rc, doc, _ = run_bash(tmp_path, body)
    assert rc == 0 and (tmp_path / "m").exists()
    assert doc["versions"]["gone"] == "not found" and doc["versions"]["bad"] == "v9.9" and doc["status"] == "ok"


def test_finalize_versions_tsv_without_tools_entry_goes_to_versions_only(tmp_path):
    tools = tmp_path / "t.tsv"; tools.write_text("flye\t0\t0\t1\n")
    vers = tmp_path / "v.tsv"; vers.write_text("flye\t2.9.6\nminimap2\t2.28\n")
    js = tmp_path / "s.json"
    r = subprocess.run([sys.executable, os.path.join(BIN, "ff_status.py"), "finalize", "--stage", "d", "--sample", "S",
                        "--json", str(js), "--tools", str(tools), "--versions", str(vers)], capture_output=True, text=True)
    doc = json.load(open(js))
    assert r.returncode == 0 and doc["tools"]["flye"]["version"] == "2.9.6" and doc["versions"]["minimap2"] == "2.28"


def test_a_stage_that_only_skipped_reports_skipped():
    """Nothing ran, so the stage did no work — `ok` would claim success for an empty output."""
    status, note = ff_status.compute_status({}, skips={"interproscan": "data not mounted"})
    assert status == "skipped"
    assert "interproscan" in note and "data not mounted" in note


def test_a_skipped_extra_alongside_real_work_is_still_ok():
    tools = {"funannotate_predict": {"exit": 0, "optional": False, "seconds": 9}}
    assert ff_status.compute_status(tools, skips={"genemark": "no licence key"}) == ("ok", "")
