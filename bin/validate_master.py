#!/usr/bin/env python3
"""Validate a master table against the master-table schema (W3.3).

  validate_master.py --master master_fungi.tsv [--schema fungiforge/resources/master_schema.json]
                     [--json report.json] [--strict]

Checks: the header carries every schema column in order (extra trailing columns are allowed and
reported, so a newer pipeline's table still validates against an older schema); `sample` is present
and unique; every value is either the NA marker or parses as its declared type; enum fields hold a
declared value. Exit 1 on any error (with --strict, also on warnings)."""
from __future__ import annotations
import argparse, csv, json, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_SCHEMA = os.path.join(HERE, "fungiforge", "resources", "master_schema.json")


def load_schema(path=None):
    p = path or DEFAULT_SCHEMA
    if not os.path.exists(p):
        try:
            from importlib.resources import files as _pkg_files
            p = str(_pkg_files("fungiforge").joinpath("resources", "master_schema.json"))
        except Exception:  # noqa: BLE001
            pass
    with open(p) as fh:
        return json.load(fh)


def check_value(value, col, na="NA"):
    """-> error string or None."""
    if value is None or value == "" or value == na:
        return None
    t = col.get("type")
    if t == "integer":
        try:
            int(float(value))
        except ValueError:
            return f"{col['name']}: '{value}' is not an integer"
        if float(value) != int(float(value)):
            return f"{col['name']}: '{value}' is not an integer"
    elif t == "number":
        try:
            float(value)
        except ValueError:
            return f"{col['name']}: '{value}' is not a number"
    elif t == "boolean":
        if str(value).lower() not in ("true", "false"):
            return f"{col['name']}: '{value}' is not a boolean"
    if col.get("enum") and value not in col["enum"]:
        return f"{col['name']}: '{value}' is not one of {col['enum']}"
    return None


def validate(master_path, schema):
    errors, warnings = [], []
    na = schema.get("na_value", "NA")
    cols = schema["columns"]
    names = [c["name"] for c in cols]
    with open(master_path) as fh:
        rows = list(csv.DictReader(fh, delimiter=schema.get("delimiter", "\t")))
        header = list(rows[0].keys()) if rows else []
    if not rows:
        with open(master_path) as fh:
            header = (fh.readline().rstrip("\n").split(schema.get("delimiter", "\t")))
        warnings.append("no data rows")
    missing = [n for n in names if n not in header]
    if missing:
        errors.append(f"missing columns: {', '.join(missing)}")
    common = [h for h in header if h in names]
    if common != [n for n in names if n in header]:
        errors.append("schema columns are out of order in the header (the master table is append-only and ordered)")
    extra = [h for h in header if h not in names]
    if extra:
        warnings.append(f"columns not in this schema version (newer pipeline?): {', '.join(extra)}")
    seen = set()
    for i, r in enumerate(rows, 1):
        s = (r.get("sample") or "").strip()
        if not s or s == na:
            errors.append(f"row {i}: empty sample")
        elif s in seen:
            errors.append(f"row {i}: duplicate sample '{s}'")
        seen.add(s)
        for c in cols:
            if c["name"] in r:
                e = check_value((r[c["name"]] or "").strip(), c, na)
                if e:
                    errors.append(f"row {i} ({s}): {e}")
    return {"master": os.path.abspath(master_path), "schema_version": schema.get("schema_version"), "n_rows": len(rows),
            "n_columns": len(header), "n_schema_columns": len(names), "errors": errors, "warnings": warnings,
            "valid": not errors}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", required=True); ap.add_argument("--schema", default=None)
    ap.add_argument("--json", default=None); ap.add_argument("--strict", action="store_true")
    a = ap.parse_args()
    rep = validate(a.master, load_schema(a.schema))
    if a.json:
        json.dump(rep, open(a.json, "w"), indent=2)
    for w in rep["warnings"]:
        print(f"[validate_master] WARNING {w}")
    for e in rep["errors"][:50]:
        print(f"[validate_master] ERROR {e}", file=sys.stderr)
    print(f"[validate_master] {rep['n_rows']} rows x {rep['n_columns']} columns against schema {rep['schema_version']}: "
          f"{'valid' if rep['valid'] else str(len(rep['errors'])) + ' error(s)'}")
    sys.exit(1 if rep["errors"] or (a.strict and rep["warnings"]) else 0)


if __name__ == "__main__":
    main()
