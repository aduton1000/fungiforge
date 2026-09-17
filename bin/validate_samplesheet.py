#!/usr/bin/env python3
"""Validate a fungiforge samplesheet before a run (W4.1, L21).

  validate_samplesheet.py --samplesheet samples.csv [--json report.json] [--strict] [--no-check-files]

Checks, with one clear message per problem (all problems at once, not the first):
  * required header `sample`, and at least one of the read columns
  * sample ids: present, unique, and safe for file names and tool arguments
  * reads: ONT alone (long-read), or BOTH Illumina mates (short-read), or both (hybrid);
    a single Illumina mate is an error; files must exist (relative paths resolve against the
    launch directory and then the samplesheet's own directory) and be non-empty
  * metadata: compartment / facility / season are free text but are reported as a distinct-value
    summary so a typo ("AIR" vs "Air") is visible before a batch runs; unknown columns are listed
Exit 1 on any error (with --strict, also on warnings)."""
from __future__ import annotations
import argparse, csv, json, os, re, sys

REQUIRED = ["sample"]
READ_COLS = ["ont_fastq", "illumina_r1", "illumina_r2"]
META_COLS = ["compartment", "facility", "season"]
SAFE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
FASTQ_SUFFIX = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fastq.bz2", ".fq.bz2")


def resolve(path, sheet_dir):
    """Absolute path, else relative to the launch dir, else next to the samplesheet."""
    if not path:
        return None
    if os.path.isabs(path):
        return path
    if os.path.exists(path):
        return os.path.abspath(path)
    near = os.path.join(sheet_dir, path)
    return near if os.path.exists(near) else os.path.abspath(path)


def validate(sheet, check_files=True, allow_pod5=True):
    errors, warnings, rows_out = [], [], []
    sheet_dir = os.path.dirname(os.path.abspath(sheet))
    with open(sheet, newline="") as fh:
        sample = fh.read(4096); fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        rdr = csv.DictReader(fh, dialect=dialect)
        header = [h.strip() for h in (rdr.fieldnames or [])]
        rows = [{(k or "").strip(): (v or "").strip() for k, v in r.items()} for r in rdr]
    for c in REQUIRED:
        if c not in header:
            errors.append(f"header: required column '{c}' is missing (found: {', '.join(header) or 'nothing'})")
    if not any(c in header for c in READ_COLS):
        errors.append(f"header: none of the read columns {READ_COLS} is present")
    unknown = [h for h in header if h not in REQUIRED + READ_COLS + META_COLS]
    if unknown:
        warnings.append(f"header: columns fungiforge ignores: {', '.join(unknown)}")
    if not rows:
        errors.append("no data rows")
    seen = {}
    for i, r in enumerate(rows, start=2):          # line 1 is the header
        sid = r.get("sample", "")
        if not sid:
            errors.append(f"line {i}: empty sample id"); continue
        if not SAFE.match(sid):
            errors.append(f"line {i}: sample id '{sid}' has characters that break file names or tool arguments "
                          "(use letters, digits, dot, dash, underscore; start with a letter or digit)")
        if sid in seen:
            errors.append(f"line {i}: duplicate sample id '{sid}' (first seen on line {seen[sid]})")
        seen[sid] = i
        ont, r1, r2 = (r.get(c, "") for c in READ_COLS)
        if bool(r1) != bool(r2):
            errors.append(f"line {i} ({sid}): Illumina must be paired — got r1='{r1}', r2='{r2}'")
        if not ont and not (r1 and r2):
            errors.append(f"line {i} ({sid}): needs ont_fastq, or both illumina_r1 and illumina_r2")
        if r1 and r2 and r1 == r2:
            errors.append(f"line {i} ({sid}): illumina_r1 and illumina_r2 are the same file")
        mode = "longread" if ont else "shortread"
        if ont and r1 and r2:
            mode = "hybrid"
        for col, raw in zip(READ_COLS, (ont, r1, r2)):
            if not raw:
                continue
            p = resolve(raw, sheet_dir)
            if check_files and not os.path.exists(p):
                errors.append(f"line {i} ({sid}): {col} not found: '{raw}' (looked in the launch directory and next to the samplesheet)")
            elif check_files and os.path.isdir(p):
                if col != "ont_fastq" or not allow_pod5:
                    errors.append(f"line {i} ({sid}): {col} is a directory: '{raw}' (only ont_fastq may be a pod5 directory, with --basecall)")
            elif check_files and os.path.getsize(p) == 0:
                errors.append(f"line {i} ({sid}): {col} is empty: '{raw}'")
            elif not p.lower().endswith(FASTQ_SUFFIX) and not (col == "ont_fastq" and os.path.isdir(p)):
                warnings.append(f"line {i} ({sid}): {col} does not look like a FASTQ ('{os.path.basename(raw)}')")
        rows_out.append({"sample": sid, "assembly_mode": mode, **{c: r.get(c, "") for c in META_COLS}})
    meta_summary = {c: sorted({(r.get(c) or "") for r in rows_out if r.get(c)}) for c in META_COLS}
    for c, vals in meta_summary.items():
        lowered = {}
        for v in vals:
            lowered.setdefault(v.strip().lower(), []).append(v)
        clashes = {k: v for k, v in lowered.items() if len(v) > 1}
        if clashes:
            warnings.append(f"{c}: values differing only in case or spacing: " +
                            "; ".join("/".join(v) for v in clashes.values()))
        if not vals and rows_out:
            warnings.append(f"{c}: empty for every isolate (the comparative layer groups by it)")
    modes = {}
    for r in rows_out:
        modes[r["assembly_mode"]] = modes.get(r["assembly_mode"], 0) + 1
    return {"samplesheet": os.path.abspath(sheet), "n_samples": len(rows_out), "modes": modes,
            "metadata_values": meta_summary, "errors": errors, "warnings": warnings, "valid": not errors}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True); ap.add_argument("--json", default=None)
    ap.add_argument("--strict", action="store_true"); ap.add_argument("--no-check-files", action="store_true")
    a = ap.parse_args()
    if not os.path.exists(a.samplesheet):
        print(f"[validate_samplesheet] ERROR samplesheet not found: {a.samplesheet}", file=sys.stderr); sys.exit(1)
    rep = validate(a.samplesheet, check_files=not a.no_check_files)
    if a.json:
        json.dump(rep, open(a.json, "w"), indent=2)
    for w in rep["warnings"]:
        print(f"[validate_samplesheet] WARNING {w}")
    for e in rep["errors"]:
        print(f"[validate_samplesheet] ERROR {e}", file=sys.stderr)
    modes = ", ".join(f"{v} {k}" for k, v in sorted(rep["modes"].items())) or "none"
    print(f"[validate_samplesheet] {rep['n_samples']} isolates ({modes}): "
          f"{'valid' if rep['valid'] else str(len(rep['errors'])) + ' error(s)'}"
          f"{'' if not rep['warnings'] else ', ' + str(len(rep['warnings'])) + ' warning(s)'}")
    sys.exit(1 if rep["errors"] or (a.strict and rep["warnings"]) else 0)


if __name__ == "__main__":
    main()
