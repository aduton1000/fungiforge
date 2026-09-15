#!/usr/bin/env python3
"""Build results/pipeline_info/provenance.json — everything needed to reproduce a run.

Two steps, because they run in different places:

  aggregate   (PROVENANCE process, inside the fungiforge image, end of the DAG)
      --run-info run_info.json      pipeline/nextflow/run metadata + effective params +
                                    process->container map (written by main.nf)
      --db-manifest db_manifest.json  output of bin/check_databases.py (DB_CHECK process)
      --jsons <every per-stage JSON>  status, exit codes, tool versions per sample/stage
      --out provenance.json

  images      (workflow.onComplete in main.nf, on the head node — where the image cache and
               the docker daemon are visible)
      --provenance provenance.json --engine apptainer|singularity|docker|none
      [--cache-dir DIR] [--work-dir DIR] [--no-hash]
      resolves every container reference to the file or image that actually ran and records
      its sha256 (files) or repo digest / image id (docker). Rewrites the file in place.

  validate    provenance.json [--schema fungiforge/resources/provenance.schema.json]
      checks the document against the shipped schema (small built-in validator: type,
      properties, required, additionalProperties, items, enum); exit 1 on violation.
"""
from __future__ import annotations
import argparse, datetime, glob, hashlib, json, os, subprocess, sys

SCHEMA_VERSION = "1.0"
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SCHEMA = os.path.join(os.path.dirname(HERE), "fungiforge", "resources", "provenance.schema.json")


# ----------------------------------------------------------------------------- aggregate
def load_json(path, what):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception as e:  # noqa: BLE001
        return {"_error": f"{what} unreadable: {e}"}


def aggregate(a):
    run_info = load_json(a.run_info, "run_info") if a.run_info else {}
    db = load_json(a.db_manifest, "db_manifest") if a.db_manifest else {}
    samples, tool_versions, conflicts, unreadable = {}, {}, {}, []
    stage_counts = {}
    paths = []
    for pat in a.jsons or []:
        paths.extend(sorted(glob.glob(pat)) if any(c in pat for c in "*?[") else [pat])
    for p in paths:
        try:
            with open(p) as fh:
                doc = json.load(fh)
        except Exception as e:  # noqa: BLE001
            unreadable.append({"file": os.path.basename(p), "error": str(e)})
            continue
        sid, stage = str(doc.get("sample", "?")), str(doc.get("stage", os.path.basename(p)))
        entry = {"status": doc.get("status", "unknown"),
                 "tools": doc.get("tools", {}),
                 "skipped_tools": doc.get("skipped_tools", {}),
                 "versions": doc.get("versions", {}),
                 "note": doc.get("note", "")}
        samples.setdefault(sid, {})[stage] = entry
        c = stage_counts.setdefault(stage, {"ok": 0, "partial": 0, "failed": 0, "skipped": 0, "unknown": 0})
        c[entry["status"] if entry["status"] in c else "unknown"] += 1
        for label, v in entry["versions"].items():
            if v in ("not found", "unknown"):
                continue
            if label in tool_versions and tool_versions[label] != v:
                conflicts.setdefault(label, sorted({tool_versions[label]})).append(v)
                conflicts[label] = sorted(set(conflicts[label]))
            else:
                tool_versions.setdefault(label, v)
    containers = run_info.get("containers", {}) if isinstance(run_info, dict) else {}
    prov = {
        "schema_version": SCHEMA_VERSION,
        "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "pipeline": run_info.get("pipeline", {}),
        "nextflow": run_info.get("nextflow", {}),
        "run": run_info.get("run", {}),
        "params": run_info.get("params", {}),
        "containers": containers if isinstance(containers, dict) else {"all": str(containers)},
        "images": {},                     # filled by `images` on the head node
        "databases": db,
        "samples": samples,
        "n_samples": len(samples),
        "stage_summary": stage_counts,
        "tool_versions": dict(sorted(tool_versions.items())),
        "version_conflicts": conflicts,   # same label, different versions across tasks — should be empty
        "unreadable_stage_json": unreadable,
    }
    with open(a.out, "w") as fh:
        json.dump(prov, fh, indent=2, sort_keys=False)
    print(f"[provenance] {len(samples)} sample(s), {len(paths)} stage JSON(s), "
          f"{len(tool_versions)} tool version(s), {len(conflicts)} conflict(s) -> {a.out}")
    if conflicts:
        print("[provenance] WARNING version conflicts: " + json.dumps(conflicts), file=sys.stderr)


# ----------------------------------------------------------------------------- images
def simple_name(ref):
    """Nextflow's cache file name for a container reference (SingularityCache.simpleName)."""
    name = ref.split("://", 1)[1] if "://" in ref else ref
    if name.endswith(".sif"):
        return name
    return name.replace(":", "-").replace("/", "-") + ".img"


def sha256_file(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def describe_file(path, do_hash):
    st = os.stat(path)
    d = {"path": path, "size": st.st_size,
         "mtime": datetime.datetime.fromtimestamp(st.st_mtime, datetime.timezone.utc).isoformat(timespec="seconds")}
    if do_hash:
        d["sha256"] = sha256_file(path)
    return d


def docker_inspect(ref):
    try:
        out = subprocess.check_output(["docker", "image", "inspect", "--format", "{{json .}}", ref],
                                      stderr=subprocess.DEVNULL, text=True)
        d = json.loads(out)
        return {"image_id": d.get("Id"), "repo_digests": d.get("RepoDigests", []),
                "created": d.get("Created"), "architecture": d.get("Architecture")}
    except Exception as e:  # noqa: BLE001
        return {"error": f"docker inspect failed: {e}"}


def resolve_image(ref, engine, cache_dirs, do_hash):
    rec = {"reference": ref, "engine": engine, "resolved": False}
    try:
        if ref.startswith("/"):                             # site .sif / .img file
            if os.path.isfile(ref):
                rec.update(kind="file", resolved=True, **describe_file(ref, do_hash))
            else:
                rec.update(kind="file", error="file not found")
        elif engine == "docker":
            info = docker_inspect(ref)
            rec.update(kind="docker", **info)
            rec["resolved"] = "error" not in info
        elif engine in ("apptainer", "singularity"):
            name = simple_name(ref)
            for d in cache_dirs:
                cand = os.path.join(d, name)
                if os.path.isfile(cand):
                    rec.update(kind="cache", resolved=True, **describe_file(cand, do_hash))
                    break
            else:
                rec["error"] = f"{name} not found in " + (", ".join(cache_dirs) or "(no cache dir known)")
        else:
            rec.update(kind="none", resolved=True, note="no container engine (conda/native run); versions come from tool_versions")
    except Exception as e:  # noqa: BLE001
        rec["error"] = str(e)
    return rec


def images(a):
    prov = load_json(a.provenance, "provenance")
    if "_error" in prov:
        print("[provenance] " + prov["_error"], file=sys.stderr)
        sys.exit(1)
    cache_dirs = [d for d in [a.cache_dir] if d]
    if a.work_dir:
        cache_dirs += [os.path.join(a.work_dir, "singularity"), os.path.join(a.work_dir, "apptainer")]
    refs = sorted({v for v in prov.get("containers", {}).values() if isinstance(v, str) and v})
    out = {}
    for ref in refs:
        out[ref] = resolve_image(ref, a.engine, cache_dirs, not a.no_hash)
    prov["images"] = out
    prov["images_hashed"] = not a.no_hash
    with open(a.out or a.provenance, "w") as fh:
        json.dump(prov, fh, indent=2)
    n_ok = sum(1 for r in out.values() if r.get("resolved"))
    print(f"[provenance] images: {n_ok}/{len(out)} resolved ({'sha256' if not a.no_hash else 'size only'}) -> {a.out or a.provenance}")
    for ref, r in out.items():
        if not r.get("resolved") and r.get("kind") != "none":
            print(f"[provenance]   unresolved {ref}: {r.get('error')}", file=sys.stderr)


# ----------------------------------------------------------------------------- validate
_TYPES = {"object": dict, "array": list, "string": str, "integer": int, "number": (int, float), "boolean": bool, "null": type(None)}


def _check(doc, schema, path, errors):
    t = schema.get("type")
    if t:
        types = t if isinstance(t, list) else [t]
        ok = any(isinstance(doc, _TYPES[x]) and not (x in ("integer", "number") and isinstance(doc, bool)) for x in types)
        if not ok:
            errors.append(f"{path or '$'}: expected {t}, got {type(doc).__name__}")
            return
    if "enum" in schema and doc not in schema["enum"]:
        errors.append(f"{path or '$'}: {doc!r} not in {schema['enum']}")
    if isinstance(doc, dict):
        for k in schema.get("required", []):
            if k not in doc:
                errors.append(f"{path or '$'}: missing required key '{k}'")
        props = schema.get("properties", {})
        for k, v in doc.items():
            if k in props:
                _check(v, props[k], f"{path}.{k}" if path else k, errors)
            elif "additionalProperties" in schema:
                ap = schema["additionalProperties"]
                if ap is False:
                    errors.append(f"{path or '$'}: unexpected key '{k}'")
                elif isinstance(ap, dict):
                    _check(v, ap, f"{path}.{k}" if path else k, errors)
    if isinstance(doc, list) and "items" in schema:
        for i, v in enumerate(doc):
            _check(v, schema["items"], f"{path}[{i}]", errors)


def validate_doc(doc, schema):
    errors = []
    _check(doc, schema, "", errors)
    return errors


def validate(a):
    doc = load_json(a.provenance, "provenance")
    schema = load_json(a.schema, "schema")
    errors = validate_doc(doc, schema)
    for e in errors:
        print("[provenance] schema violation: " + e, file=sys.stderr)
    print(f"[provenance] {a.provenance}: {'valid' if not errors else str(len(errors)) + ' violation(s)'}")
    sys.exit(1 if errors else 0)


# ----------------------------------------------------------------------------- cli
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("aggregate"); g.add_argument("--run-info"); g.add_argument("--db-manifest")
    g.add_argument("--jsons", nargs="*", default=[]); g.add_argument("--out", default="provenance.json"); g.set_defaults(fn=aggregate)
    i = sub.add_parser("images"); i.add_argument("--provenance", required=True); i.add_argument("--engine", default="none")
    i.add_argument("--cache-dir", default=""); i.add_argument("--work-dir", default=""); i.add_argument("--no-hash", action="store_true")
    i.add_argument("--out"); i.set_defaults(fn=images)
    v = sub.add_parser("validate"); v.add_argument("provenance"); v.add_argument("--schema", default=DEFAULT_SCHEMA); v.set_defaults(fn=validate)
    a = ap.parse_args(); a.fn(a)


if __name__ == "__main__":
    main()
