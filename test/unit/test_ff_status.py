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
    rc, doc, _ = run_bash(tmp_path, "ff_skip genomad 'database not staged'")
    assert rc == 0 and doc["status"] == "ok" and doc["skipped_tools"]["genomad"] == "database not staged"
