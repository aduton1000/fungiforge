#!/usr/bin/env python3
"""Verify the reference-database root before a run starts (Stage 00: DB_CHECK).

Each database fungiforge uses has a completion marker (`.done`, written by
bin/fetch_references.sh) and at least one file that must exist. Missing or incomplete
databases are reported together, once, before any compute is spent; the run fails unless
--allow-missing is given (in which case dependent stages will record `skipped`).

Also writes a manifest (JSON) of what was found — name, marker timestamp, key files with
sizes, MANIFEST.tsv source line — which the PROVENANCE stage embeds in provenance.json.

  check_databases.py --data-dir DIR --out db_manifest.json [--require a,b,c] [--allow-missing]
                     [--path funannotate=/other/dir ...]   # honour --funannotate_db / --antismash_db overrides
"""
from __future__ import annotations
import argparse, glob, json, os, sys

# name -> (required key file globs relative to <data_dir>/<name>, human description)
DATABASES = {
    "unite":        (["sh_general_release*dynamic*.fasta"], "UNITE ITS reference (stage 08)"),
    "kraken2":      (["hash.k2d", "taxo.k2d", "opts.k2d"], "Kraken2 PlusPF (stage 04, triage)"),
    "busco":        (["lineages/fungi_odb10", "fungi_odb10", "busco_downloads/lineages/fungi_odb10"], "BUSCO fungi_odb10 (stage 05; any one path)"),
    "funannotate":  (["Pfam-A.hmm", "funannotate-db-info.txt"], "funannotate databases (stage 07)"),
    "antismash":    (["clusterblast", "knownclusterblast", "pfam"], "antiSMASH databases (stage 11)"),
    "fungamr":      (["FungAMR_*.tsv", "reference_proteins.faa"], "FungAMR tables + reference proteins (stage 09)"),
    "refseq_fungi": (["*.zip"], "sourmash fungal genome signatures (stage 08, --genome_id)"),
}
# Present-or-not only; stages record `skipped` for these when absent.
OPTIONAL = {
    "eggnog":               "eggNOG (stage 07, used when present)",
    "genomad_db":           "geNomad (stage 10, used when present)",
    "refseq_fungi_genomes": "reference genome set for ANI novelty (stage 12, used when present)",
    "rvdb":                 "RVDB-prot mycovirus references (fetched; consumed from W2.8 onwards)",
}


def read_manifest(data_dir):
    rows = {}
    p = os.path.join(data_dir, "MANIFEST.tsv")
    if os.path.exists(p):
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0] != "database":
                rows[f[0]] = {"detail": f[1], "status": f[2], "timestamp": f[3]}   # last line wins
    return rows


def check_one(data_dir, name, key_globs, any_of=False, override=None):
    d = override or os.path.join(data_dir, name)
    res = {"path": d, "present": os.path.isdir(d), "done_marker": None, "key_files": {}, "missing": []}
    if not res["present"]:
        res["missing"] = ["directory"]
        return res
    marker = os.path.join(d, ".done")
    res["done_marker"] = open(marker).read().strip() if os.path.exists(marker) else None
    hits = 0
    for g in key_globs:
        found = glob.glob(os.path.join(d, g))
        if found:
            hits += 1
            for f in found[:5]:
                res["key_files"][os.path.relpath(f, d)] = os.path.getsize(f) if os.path.isfile(f) else "dir"
        else:
            res["missing"].append(g)
    if any_of and hits:            # e.g. busco: one of several layouts is enough
        res["missing"] = []
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--require", default=",".join(DATABASES), help="comma list of required DB names")
    ap.add_argument("--allow-missing", action="store_true")
    ap.add_argument("--path", action="append", default=[], metavar="NAME=DIR", help="use DIR instead of <data_dir>/NAME")
    a = ap.parse_args()
    overrides = dict(p.split("=", 1) for p in a.path if "=" in p and p.split("=", 1)[1])

    required = [x for x in a.require.split(",") if x]
    manifest_rows = read_manifest(a.data_dir)
    report = {"data_dir": a.data_dir, "manifest_tsv": os.path.exists(os.path.join(a.data_dir, "MANIFEST.tsv")),
              "databases": {}, "optional": {}, "problems": []}
    for name in required:
        globs, desc = DATABASES.get(name, ([], name))
        r = check_one(a.data_dir, name, globs, any_of=(name == "busco"), override=overrides.get(name))
        r["description"] = desc
        r["source"] = manifest_rows.get(name)
        if r["missing"]:
            report["problems"].append(f"{name}: missing {', '.join(r['missing'])} ({desc})")
        elif r["done_marker"] is None:
            report["problems"].append(f"{name}: no .done marker — fetch may be incomplete ({desc})")
        report["databases"][name] = r
    for name, desc in OPTIONAL.items():
        d = os.path.join(a.data_dir, name)
        present = os.path.isdir(d) and bool(os.listdir(d))
        report["optional"][name] = {"present": present, "description": desc, "source": manifest_rows.get(name)}
    report["ok"] = not report["problems"]
    json.dump(report, open(a.out, "w"), indent=2)

    for name, r in report["databases"].items():
        print(f"[db_check] {name:14s} {'OK     ' if not r['missing'] and r['done_marker'] else 'PROBLEM'}  {r['description']}")
    for name, r in report["optional"].items():
        print(f"[db_check] {name:14s} {'present' if r['present'] else 'absent '}  (optional) {r['description']}")
    if report["problems"]:
        print("[db_check] problems:\n  - " + "\n  - ".join(report["problems"]), file=sys.stderr)
        if not a.allow_missing:
            print("[db_check] refusing to start: fix the databases (bin/fetch_references.sh) or pass --allow_missing_db", file=sys.stderr)
            sys.exit(1)
    print(f"[db_check] manifest -> {a.out}")


if __name__ == "__main__":
    main()
