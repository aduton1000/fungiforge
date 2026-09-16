#!/usr/bin/env python3
"""Compare a finished fungiforge run against an expected-results file (W1.1 validation suite).

  validate_run.py --results results/ --expected test/expected/AfumCEA10.json [--sample ID]
                  [--out validation.tsv] [--json validation.json]

Expected file (JSON):
  {"sample": "AfumCEA10", "source": {...free text...},
   "checks": {
     "assembly_len":   {"from": "assembly_qc.assembly_len", "expected": 29900000, "rel_tol": 0.02},
     "n_contigs":      {"from": "assembly_qc.n_contigs", "max": 60},
     "busco_complete": {"from": "assembly_qc.busco_complete", "min": 98.0},
     "species":        {"from": "identify.species", "equals": "Aspergillus fumigatus"},
     "verdict":        {"from": "decontam.verdict", "one_of": ["fungal", "likely_fungal"]},
     "resistance":     {"from": "master.resistant_classes", "equals": "none"},
     "id_genus":       {"from": "identify.species", "regex": "^Aspergillus "},
     "stages":         {"from": "master.stages_failed", "equals": "none", "optional": false}
   }}

`from` names where the observed value lives:
  <stage>.<key>[.<key>...]   the stage JSON  results/<sample>/*/<sample>.<stage>.json
  master.<column>            results/04_summary/<sample>.master.tsv
  provenance.<key>[...]      results/pipeline_info/provenance.json

Rules (any combination; all must hold): equals, one_of, regex, min, max, and
expected with abs_tol and/or rel_tol. A check whose value cannot be found is MISSING and
counts as a failure unless "optional": true. Exit status 1 when any check fails.
"""
from __future__ import annotations
import argparse, csv, glob, json, os, re, sys


def load_stage_json(results, sample, stage):
    hits = glob.glob(os.path.join(results, sample, "*", f"{sample}.{stage}.json"))
    if not hits:
        return None
    with open(sorted(hits)[0]) as fh:
        return json.load(fh)


def load_master(results, sample):
    p = os.path.join(results, "04_summary", f"{sample}.master.tsv")
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else None


def load_provenance(results):
    p = os.path.join(results, "pipeline_info", "provenance.json")
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


def dig(doc, keys):
    cur = doc
    for k in keys:
        if isinstance(cur, dict) and k in cur:
            cur = cur[k]
        elif isinstance(cur, list) and k.isdigit() and int(k) < len(cur):
            cur = cur[int(k)]
        else:
            return None, False
    return cur, True


class Run:
    """Lazy access to a run's documents by `from` path."""

    def __init__(self, results, sample):
        self.results, self.sample, self._cache = results, sample, {}

    def get(self, frm):
        parts = frm.split(".")
        src, keys = parts[0], parts[1:]
        if src not in self._cache:
            if src == "master":
                self._cache[src] = load_master(self.results, self.sample)
            elif src == "provenance":
                self._cache[src] = load_provenance(self.results)
            else:
                self._cache[src] = load_stage_json(self.results, self.sample, src)
        doc = self._cache[src]
        if doc is None:
            return None, False
        return dig(doc, keys)


def coerce_number(v):
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        try:
            return float(v)
        except ValueError:
            return None
    return None


def evaluate(rule, observed, found):
    """Return (status, detail) for one check. status in PASS | FAIL | MISSING."""
    if not found or observed in (None, "", "NA"):
        return ("PASS" if rule.get("optional") else "MISSING"), "value not found"
    problems = []
    if "equals" in rule:
        exp = rule["equals"]
        same = (str(observed).lower() == str(exp).lower()) if isinstance(exp, (str, bool)) else (coerce_number(observed) == coerce_number(exp))
        if not same:
            problems.append(f"expected == {exp!r}")
    if "one_of" in rule and str(observed) not in [str(x) for x in rule["one_of"]]:
        problems.append(f"expected one of {rule['one_of']}")
    if "regex" in rule and not re.search(rule["regex"], str(observed)):
        problems.append(f"expected to match /{rule['regex']}/")
    num = coerce_number(observed)
    for key in ("min", "max", "expected"):
        if key in rule and num is None:
            problems.append(f"{key} rule needs a number, got {observed!r}")
    if num is not None:
        if "min" in rule and num < float(rule["min"]):
            problems.append(f"expected >= {rule['min']}")
        if "max" in rule and num > float(rule["max"]):
            problems.append(f"expected <= {rule['max']}")
        if "expected" in rule and coerce_number(rule["expected"]) is not None:
            exp = float(rule["expected"])
            tol = max(float(rule.get("abs_tol", 0)), abs(exp) * float(rule.get("rel_tol", 0)))
            if abs(num - exp) > tol:
                problems.append(f"expected {exp:g} ± {tol:g}")
    return ("FAIL" if problems else "PASS"), "; ".join(problems)


def describe_rule(rule):
    bits = []
    for k in ("equals", "one_of", "regex", "min", "max"):
        if k in rule:
            bits.append(f"{k}={rule[k]!r}")
    if "expected" in rule:
        t = []
        if "abs_tol" in rule: t.append(f"±{rule['abs_tol']}")
        if "rel_tol" in rule: t.append(f"±{100*float(rule['rel_tol']):g}%")
        bits.append(f"expected={rule['expected']} {' '.join(t)}".strip())
    if rule.get("optional"):
        bits.append("optional")
    return " ".join(bits)


def validate(results, expected, sample=None):
    sample = sample or expected["sample"]
    run = Run(results, sample)
    rows = []
    for name, rule in expected["checks"].items():
        observed, found = run.get(rule["from"])
        status, detail = evaluate(rule, observed, found)
        rows.append({"check": name, "from": rule["from"], "observed": observed if found else None,
                     "rule": describe_rule(rule), "status": status, "detail": detail})
    n_fail = sum(r["status"] in ("FAIL", "MISSING") for r in rows)
    return {"sample": sample, "results": os.path.abspath(results), "expected_source": expected.get("source", {}),
            "n_checks": len(rows), "n_fail": n_fail, "passed": n_fail == 0, "checks": rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--expected", required=True)
    ap.add_argument("--sample", default=None, help="sample id (default: the expected file's)")
    ap.add_argument("--out", default=None, help="write the table as TSV")
    ap.add_argument("--json", default=None, help="write the full result as JSON")
    a = ap.parse_args()
    with open(a.expected) as fh:
        expected = json.load(fh)
    rep = validate(a.results, expected, a.sample)
    w = max(len(r["check"]) for r in rep["checks"]) if rep["checks"] else 5
    print(f"validation: {rep['sample']}  ({rep['results']})")
    for r in rep["checks"]:
        obs = r["observed"] if r["observed"] is not None else "-"
        print(f"  {r['status']:7s} {r['check']:{w}s}  observed={obs!s:<28.28} {r['rule']}" + (f"  [{r['detail']}]" if r["detail"] and r["status"] != "PASS" else ""))
    print(f"  {'PASS' if rep['passed'] else 'FAIL'}: {rep['n_checks'] - rep['n_fail']}/{rep['n_checks']} checks passed")
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("sample\tcheck\tfrom\tobserved\trule\tstatus\tdetail\n")
            for r in rep["checks"]:
                fh.write("\t".join(str(x) for x in (rep["sample"], r["check"], r["from"], r["observed"], r["rule"], r["status"], r["detail"])) + "\n")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=2)
    sys.exit(0 if rep["passed"] else 1)


if __name__ == "__main__":
    main()
