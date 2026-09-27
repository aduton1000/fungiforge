"""Unit tests for bin/make_provenance.py (aggregate / images / validate) and the shipped schema."""
import hashlib, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
BIN = os.path.join(ROOT, "bin")
SCRIPT = os.path.join(BIN, "make_provenance.py")
SCHEMA = os.path.join(ROOT, "fungiforge", "resources", "provenance.schema.json")
sys.path.insert(0, BIN)
import make_provenance  # noqa: E402

RUN_INFO = {
    "pipeline": {"name": "fungiforge", "version": "0.2.0", "repository": None, "revision": None,
                 "commit_id": "abc123", "git_dirty": False, "script_id": "deadbeef", "project_dir": "/repo"},
    "nextflow": {"version": "26.04.3", "build": 12259},
    "run": {"session_id": "s", "run_name": "r", "start": "2026-09-15T00:00:00", "command_line": "nextflow run main.nf",
            "launch_dir": "/l", "work_dir": "/w", "profile": "local,test", "container_engine": "apptainer",
            "stub_run": False, "resume": False, "user": "u", "host": "h"},
    "containers": {"withLabel:assembly": "staphb/flye:2.9.6", "withLabel:report": "/site/images/fungiforge-0.2.0.sif"},
    "params": {"outdir": "results", "data_dir": "/db"},
}


def stage_json(tmp_path, sample, stage, status="ok", versions=None, tools=None):
    p = tmp_path / f"{sample}.{stage}.json"
    p.write_text(json.dumps({"sample": sample, "stage": stage, "status": status,
                             "tools": tools or {}, "skipped_tools": {}, "versions": versions or {}}))
    return p


def run_cli(*args):
    return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)


def aggregate(tmp_path, jsons, db=None):
    ri = tmp_path / "run_info.json"; ri.write_text(json.dumps(RUN_INFO))
    dbm = tmp_path / "db.json"; dbm.write_text(json.dumps(db or {"ok": True, "databases": {}}))
    out = tmp_path / "provenance.json"
    r = run_cli("aggregate", "--run-info", str(ri), "--db-manifest", str(dbm), "--jsons", *map(str, jsons), "--out", str(out))
    assert r.returncode == 0, r.stderr
    return json.load(open(out)), r


# ---------------------------------------------------------------- aggregate
def test_aggregate_collects_samples_stages_and_versions(tmp_path):
    js = [stage_json(tmp_path, "S1", "assemble", versions={"flye": "2.9.6-b1802"}, tools={"flye": {"exit": 0, "optional": False, "seconds": 3, "version": "2.9.6-b1802"}}),
          stage_json(tmp_path, "S1", "medaka", status="partial", versions={"medaka": "medaka 2.2.2"}),
          stage_json(tmp_path, "S2", "assemble", versions={"flye": "2.9.6-b1802"})]
    doc, r = aggregate(tmp_path, js)
    assert doc["n_samples"] == 2 and set(doc["samples"]) == {"S1", "S2"}
    assert doc["samples"]["S1"]["medaka"]["status"] == "partial"
    assert doc["stage_summary"]["assemble"] == {"ok": 2, "partial": 0, "failed": 0, "skipped": 0, "unknown": 0}
    assert doc["tool_versions"] == {"flye": "2.9.6-b1802", "medaka": "medaka 2.2.2"}
    assert doc["version_conflicts"] == {} and doc["unreadable_stage_json"] == []
    assert doc["pipeline"]["commit_id"] == "abc123" and doc["containers"]["withLabel:assembly"] == "staphb/flye:2.9.6"
    assert doc["params"]["data_dir"] == "/db" and doc["databases"]["ok"] is True


def test_aggregate_flags_version_conflicts_and_ignores_not_found(tmp_path):
    js = [stage_json(tmp_path, "S1", "assemble", versions={"flye": "2.9.6", "minimap2": "not found"}),
          stage_json(tmp_path, "S2", "assemble", versions={"flye": "2.9.5", "minimap2": "2.28-r1209"})]
    doc, r = aggregate(tmp_path, js)
    assert doc["version_conflicts"] == {"flye": ["2.9.5", "2.9.6"]}
    assert doc["tool_versions"]["minimap2"] == "2.28-r1209"
    assert "version conflicts" in r.stderr


def test_aggregate_records_unreadable_stage_json_instead_of_crashing(tmp_path):
    bad = tmp_path / "S1.broken.json"; bad.write_text("{oops")
    doc, _ = aggregate(tmp_path, [stage_json(tmp_path, "S1", "assemble"), bad])
    assert doc["unreadable_stage_json"] == [{"file": "S1.broken.json", "error": doc["unreadable_stage_json"][0]["error"]}]
    assert doc["n_samples"] == 1


def test_aggregate_expands_globs(tmp_path):
    stage_json(tmp_path, "S1", "assemble"); stage_json(tmp_path, "S1", "medaka")
    doc, _ = aggregate(tmp_path, [str(tmp_path / "S1.*.json")])
    assert set(doc["samples"]["S1"]) == {"assemble", "medaka"}


# ---------------------------------------------------------------- images
def test_simple_name_matches_nextflow_cache_naming():
    assert make_provenance.simple_name("staphb/flye:2.9.6") == "staphb-flye-2.9.6.img"
    assert make_provenance.simple_name("docker://ezlabgva/busco:v5.7.1_cv1") == "ezlabgva-busco-v5.7.1_cv1.img"
    d = "nextgenusfs/funannotate@sha256:2eadfeb46368052d8862ff3851a8af88f3468308824245fe83de68669fb3803b"
    assert make_provenance.simple_name(d) == "nextgenusfs-funannotate@sha256-2eadfeb46368052d8862ff3851a8af88f3468308824245fe83de68669fb3803b.img"
    assert make_provenance.simple_name("/site/x.sif") == "/site/x.sif"


def test_images_resolves_site_sif_and_cached_img_with_sha256(tmp_path):
    sif = tmp_path / "fungiforge-0.2.0.sif"; sif.write_bytes(b"SIF" * 1000)
    cache = tmp_path / "cache"; cache.mkdir()
    img = cache / "staphb-flye-2.9.6.img"; img.write_bytes(b"IMG" * 500)
    doc, _ = aggregate(tmp_path, [stage_json(tmp_path, "S1", "assemble")])
    prov = tmp_path / "provenance.json"
    d = json.load(open(prov)); d["containers"] = {"withLabel:assembly": "staphb/flye:2.9.6", "withLabel:report": str(sif),
                                                  "withLabel:busco": "ezlabgva/busco:v5.7.1_cv1"}
    json.dump(d, open(prov, "w"))
    r = run_cli("images", "--provenance", str(prov), "--engine", "apptainer", "--cache-dir", str(cache))
    assert r.returncode == 0, r.stderr
    d = json.load(open(prov))
    assert d["images_hashed"] is True
    assert d["images"][str(sif)]["kind"] == "file" and d["images"][str(sif)]["sha256"] == hashlib.sha256(b"SIF" * 1000).hexdigest()
    assert d["images"]["staphb/flye:2.9.6"]["kind"] == "cache" and d["images"]["staphb/flye:2.9.6"]["path"] == str(img)
    assert d["images"]["staphb/flye:2.9.6"]["sha256"] == hashlib.sha256(b"IMG" * 500).hexdigest()
    busco = d["images"]["ezlabgva/busco:v5.7.1_cv1"]
    assert busco["resolved"] is False and "not found" in busco["error"]
    assert "unresolved ezlabgva/busco" in r.stderr


def test_images_no_hash_records_size_only(tmp_path):
    sif = tmp_path / "x.sif"; sif.write_bytes(b"x" * 10)
    doc, _ = aggregate(tmp_path, [stage_json(tmp_path, "S1", "assemble")])
    prov = tmp_path / "provenance.json"
    d = json.load(open(prov)); d["containers"] = {"withLabel:report": str(sif)}; json.dump(d, open(prov, "w"))
    r = run_cli("images", "--provenance", str(prov), "--engine", "apptainer", "--no-hash")
    d = json.load(open(prov))
    assert r.returncode == 0 and d["images_hashed"] is False
    assert d["images"][str(sif)]["size"] == 10 and "sha256" not in d["images"][str(sif)]


def test_images_without_engine_marks_none(tmp_path):
    doc, _ = aggregate(tmp_path, [stage_json(tmp_path, "S1", "assemble")])
    prov = tmp_path / "provenance.json"
    r = run_cli("images", "--provenance", str(prov), "--engine", "none")
    d = json.load(open(prov))
    assert r.returncode == 0
    assert all(v["kind"] == "none" and v["resolved"] for k, v in d["images"].items() if not k.startswith("/"))
    assert d["images"]["/site/images/fungiforge-0.2.0.sif"] == {"reference": "/site/images/fungiforge-0.2.0.sif", "engine": "none",
                                                                "resolved": False, "kind": "file", "error": "file not found"}


# ---------------------------------------------------------------- validate / schema
def test_aggregated_and_imaged_document_validates_against_schema(tmp_path):
    sif = tmp_path / "x.sif"; sif.write_bytes(b"x")
    doc, _ = aggregate(tmp_path, [stage_json(tmp_path, "S1", "assemble", versions={"flye": "2.9.6"})])
    prov = tmp_path / "provenance.json"
    d = json.load(open(prov)); d["containers"]["withLabel:report"] = str(sif); json.dump(d, open(prov, "w"))
    run_cli("images", "--provenance", str(prov), "--engine", "apptainer")
    r = run_cli("validate", str(prov))
    assert r.returncode == 0 and "valid" in r.stdout, r.stderr


def test_validate_rejects_missing_keys_and_bad_status(tmp_path):
    doc, _ = aggregate(tmp_path, [stage_json(tmp_path, "S1", "assemble", status="great")])
    prov = tmp_path / "provenance.json"
    d = json.load(open(prov)); del d["tool_versions"]; json.dump(d, open(prov, "w"))
    r = run_cli("validate", str(prov))
    assert r.returncode == 1
    assert "missing required key 'tool_versions'" in r.stderr and "'great' not in" in r.stderr


def test_mini_validator_semantics():
    schema = {"type": "object", "required": ["a"], "properties": {"a": {"type": ["integer", "null"]}, "b": {"type": "array", "items": {"type": "string"}}},
              "additionalProperties": False}
    assert make_provenance.validate_doc({"a": 1, "b": ["x"]}, schema) == []
    assert make_provenance.validate_doc({"a": None}, schema) == []
    errs = make_provenance.validate_doc({"a": True, "b": [1], "c": 0}, schema)
    assert any("expected ['integer', 'null']" in e for e in errs)      # bool is not an integer
    assert any("b[0]" in e for e in errs) and any("unexpected key 'c'" in e for e in errs)



def test_images_marks_config_only_containers_as_unused_from_the_trace(tmp_path):
    sif = tmp_path / "site.sif"; sif.write_bytes(b"x" * 10)
    doc, _ = aggregate(tmp_path, [stage_json(tmp_path, "S1", "assemble")])
    prov = tmp_path / "provenance.json"
    d = json.load(open(prov)); d["containers"] = {"withLabel:report": "aduton1000/fungiforge:0.1.0", "withLabel:report|readqc": str(sif)}
    json.dump(d, open(prov, "w"))
    (tmp_path / "trace.txt").write_text("task_id\thash\tname\tstatus\tcontainer\n1\tab/cd\tREPORT (S1)\tCOMPLETED\t" + str(sif) + "\n")
    r = run_cli("images", "--provenance", str(prov), "--engine", "apptainer", "--no-hash")
    assert r.returncode == 0, r.stderr
    d = json.load(open(prov))
    assert d["images"][str(sif)]["used"] is True and d["images"][str(sif)]["resolved"] is True
    assert d["images"]["aduton1000/fungiforge:0.1.0"]["used"] is False and "overridden" in d["images"]["aduton1000/fungiforge:0.1.0"]["note"]
    assert "1/1 used images resolved" in r.stdout and "1 listed but unused" in r.stdout
