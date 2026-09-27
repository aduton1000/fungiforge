#!/usr/bin/env python3
"""Stage the type-material reference sets for the secondary identification loci (W2.3).

For each locus an NCBI Entrez query restricted to fungal records flagged "sequence from type"
(the INSDC type-material annotation) is run and the FASTA records are written to
<out_dir>/<locus>.fasta. Only records from type material are used so that a hit names the
species its type strain was described under, not whatever a submitter called an isolate.

  fetch_marker_refs.py --out-dir <data_dir>/markers [--loci CaM,BenA,TEF1,RPB2,LSU] [--batch 5000]

NCBI_API_KEY (optional) raises the request rate; NCBI_EMAIL is sent as the contact address.
Idempotent per locus: an existing non-empty <locus>.fasta is kept.
"""
from __future__ import annotations
import argparse, json, os, sys, time, urllib.parse, urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
QUERIES = {
    "CaM":  '(calmodulin[Title] OR cmdA[Title] OR CaM[Title]) AND fungi[Organism] AND "sequence from type"[Filter] NOT unverified[Title] AND 200:5000[Sequence Length]',
    "BenA": '(beta-tubulin[Title] OR "beta tubulin"[Title] OR benA[Title] OR tub2[Title]) AND fungi[Organism] AND "sequence from type"[Filter] NOT unverified[Title] AND 200:5000[Sequence Length]',
    "TEF1": '("elongation factor 1-alpha"[Title] OR "elongation factor 1 alpha"[Title] OR "translation elongation factor"[Title] OR tef1[Title] OR TEF-1[Title]) AND fungi[Organism] AND "sequence from type"[Filter] NOT unverified[Title] AND 200:5000[Sequence Length]',
    "RPB2": '("RNA polymerase II"[Title] OR rpb2[Title] OR RPB2[Title]) AND fungi[Organism] AND "sequence from type"[Filter] NOT unverified[Title] AND 200:5000[Sequence Length]',
    "LSU":  '("28S ribosomal RNA"[Title] OR "large subunit ribosomal RNA"[Title] OR "28S rRNA"[Title] OR "D1/D2"[Title]) AND fungi[Organism] AND "sequence from type"[Filter] NOT unverified[Title] AND 300:6000[Sequence Length]',
}


def _get(url, params, retries=5, sleep=1.0):
    q = dict(params)
    q["tool"] = "fungiforge"
    if os.environ.get("NCBI_EMAIL"):
        q["email"] = os.environ["NCBI_EMAIL"]
    if os.environ.get("NCBI_API_KEY"):
        q["api_key"] = os.environ["NCBI_API_KEY"]
    full = url + "?" + urllib.parse.urlencode(q)
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(full, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(sleep * (2 ** i))
    raise RuntimeError(f"NCBI request failed after {retries} tries: {full} ({last})")


def esearch(term, db="nuccore"):
    doc = json.loads(_get(f"{EUTILS}/esearch.fcgi", {"db": db, "term": term, "usehistory": "y", "retmax": 0, "retmode": "json"}))
    r = doc["esearchresult"]
    return int(r["count"]), r["webenv"], r["querykey"]


def efetch_fasta(webenv, qkey, count, batch=5000, db="nuccore", delay=0.4):
    for start in range(0, count, batch):
        yield _get(f"{EUTILS}/efetch.fcgi", {"db": db, "WebEnv": webenv, "query_key": qkey, "retstart": start,
                                             "retmax": batch, "rettype": "fasta", "retmode": "text"})
        time.sleep(delay)


def fetch_locus(locus, term, out_path, batch=5000, log=print):
    count, webenv, qkey = esearch(term)
    log(f"[markers] {locus}: {count} type-material records")
    n = 0
    tmp = out_path + ".part"
    with open(tmp, "w") as fh:
        for chunk in efetch_fasta(webenv, qkey, count, batch):
            n += chunk.count("\n>") + (1 if chunk.startswith(">") else 0)
            fh.write(chunk if chunk.endswith("\n") else chunk + "\n")
    os.replace(tmp, out_path)
    return count, n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--loci", default=",".join(QUERIES))
    ap.add_argument("--batch", type=int, default=5000)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    summary, failed = {}, []
    for locus in [x for x in a.loci.split(",") if x]:
        if locus not in QUERIES:
            print(f"[markers] unknown locus {locus}; known: {', '.join(QUERIES)}", file=sys.stderr); failed.append(locus); continue
        out = os.path.join(a.out_dir, f"{locus}.fasta")
        if os.path.exists(out) and os.path.getsize(out) > 0 and not a.force:
            print(f"[markers] {locus}: present — skip"); summary[locus] = {"records": None, "kept": True}; continue
        try:
            count, n = fetch_locus(locus, QUERIES[locus], out, a.batch)
            summary[locus] = {"records": count, "written": n, "query": QUERIES[locus]}
            print(f"[markers] {locus}: wrote {n} sequences -> {out}")
        except Exception as e:  # noqa: BLE001
            print(f"[markers] {locus}: FAILED ({e})", file=sys.stderr); failed.append(locus)
    with open(os.path.join(a.out_dir, "markers_manifest.json"), "w") as fh:
        json.dump({"date": time.strftime("%Y-%m-%d"), "loci": summary, "failed": failed}, fh, indent=2)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
