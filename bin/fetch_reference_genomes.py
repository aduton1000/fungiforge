#!/usr/bin/env python3
"""Stage the reference genome set for genome-ANI novelty (W2.9, stage 12).

For every genus of interest the NCBI Datasets API lists the reference/representative genome of
each species (`filters.reference_only`), plus any accessions given explicitly; each is downloaded
as FASTA to <out_dir>/<accession>.fna with a manifest (accession, organism, species, strain,
category, level, length, source genus). Idempotent: present files are kept.

  fetch_reference_genomes.py --out-dir <data_dir>/refseq_fungi_genomes
        [--genera Aspergillus,Candida,...] [--genera-file fungiforge/resources/novelty_genera.txt]
        [--accessions GCF_000002655.1,...] [--max-per-genus 0] [--dry-run]

NCBI_API_KEY (optional) raises the request rate. The default genus list (fungiforge/resources/
novelty_genera.txt) covers the clinically and environmentally relevant fungi; the whole fungal
reference set (~5,500 genomes) is possible with --genera Fungi but is >150 GB.
"""
from __future__ import annotations
import argparse, csv, io, json, os, re, sys, time, urllib.parse, urllib.request, zipfile

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"


def _headers():
    h = {"Accept": "application/json", "User-Agent": "fungiforge"}
    if os.environ.get("NCBI_API_KEY"):
        h["api-key"] = os.environ["NCBI_API_KEY"]
    return h


def _get(url, retries=5, binary=False):
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=_headers()), timeout=300) as r:
                data = r.read()
                return data if binary else data.decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e; time.sleep(2 ** i)
    raise RuntimeError(f"request failed: {url} ({last})")


def list_reference_genomes(taxon, max_n=0):
    """Reference/representative genome per species under `taxon` (name or taxid)."""
    out, token = [], None
    while True:
        q = {"filters.reference_only": "true", "page_size": 100}
        if token:
            q["page_token"] = token
        url = f"{API}/genome/taxon/{urllib.parse.quote(str(taxon))}/dataset_report?" + urllib.parse.urlencode(q)
        doc = json.loads(_get(url))
        for r in doc.get("reports", []):
            org = r.get("organism", {})
            name = org.get("organism_name", "")
            sp = " ".join(name.split()[:2]) if len(name.split()) >= 2 else name
            info = r.get("assembly_info", {})
            out.append({"accession": r["accession"], "organism": name, "species": sp, "strain": (org.get("infraspecific_names") or {}).get("strain", ""),
                        "tax_id": org.get("tax_id"), "category": info.get("refseq_category", ""), "level": info.get("assembly_level", ""),
                        "length": (r.get("assembly_stats") or {}).get("total_sequence_length", ""), "source": str(taxon)})
            if max_n and len(out) >= max_n:
                return out
        token = doc.get("next_page_token")
        if not token:
            break
    return out


def accession_report(accession):
    """Metadata of one explicit accession (organism, species, strain, category); blank fields when the lookup fails."""
    row = {"accession": accession, "organism": "", "species": "", "strain": "", "tax_id": "", "category": "explicit", "level": "", "length": "", "source": "accession"}
    try:
        doc = json.loads(_get(f"{API}/genome/accession/{urllib.parse.quote(accession)}/dataset_report"))
        r = (doc.get("reports") or [None])[0]
        if r:
            org = r.get("organism", {}); name = org.get("organism_name", ""); info = r.get("assembly_info", {})
            row.update(organism=name, species=" ".join(name.split()[:2]), strain=(org.get("infraspecific_names") or {}).get("strain", ""),
                       tax_id=org.get("tax_id", ""), category=info.get("refseq_category", "") or "explicit", level=info.get("assembly_level", ""),
                       length=(r.get("assembly_stats") or {}).get("total_sequence_length", ""))
    except Exception as e:  # noqa: BLE001
        print(f"[genomes] {accession}: metadata lookup failed ({e})", file=sys.stderr)
    return row


def download_fasta(accession, out_path, log=print):
    url = f"{API}/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA&filename={accession}.zip"
    data = _get(url, binary=True)
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        fnas = [n for n in z.namelist() if n.endswith(".fna")]
        if not fnas:
            raise RuntimeError(f"{accession}: no .fna in the download")
        tmp = out_path + ".part"
        with open(tmp, "wb") as fh:
            for n in fnas:
                fh.write(z.read(n))
        os.replace(tmp, out_path)
    return len(fnas)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True); ap.add_argument("--genera", default=None); ap.add_argument("--genera-file", default=None)
    ap.add_argument("--accessions", default=""); ap.add_argument("--max-per-genus", type=int, default=0); ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    genera = [g.strip() for g in (a.genera or "").split(",") if g.strip()]
    if not genera and a.genera_file and os.path.exists(a.genera_file):
        genera = [l.strip() for l in open(a.genera_file) if l.strip() and not l.startswith("#")]
    rows = []
    for g in genera:
        try:
            got = list_reference_genomes(g, a.max_per_genus)
            print(f"[genomes] {g}: {len(got)} reference genomes"); rows.extend(got)
        except Exception as e:  # noqa: BLE001
            print(f"[genomes] {g}: listing FAILED ({e})", file=sys.stderr)
    for acc in [x.strip() for x in a.accessions.split(",") if x.strip()]:
        if not any(r["accession"] == acc for r in rows):
            rows.append(accession_report(acc))
    seen, uniq = set(), []
    for r in rows:
        if r["accession"] not in seen:
            seen.add(r["accession"]); uniq.append(r)
    man = os.path.join(a.out_dir, "manifest.tsv")
    with open(man, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "organism", "species", "strain", "tax_id", "category", "level", "length", "source", "file"], delimiter="\t")
        w.writeheader()
        for r in uniq:
            w.writerow(dict(r, file=f"{r['accession']}.fna"))
    print(f"[genomes] {len(uniq)} genomes listed -> {man}")
    if a.dry_run:
        return
    failed = 0
    for i, r in enumerate(uniq, 1):
        out = os.path.join(a.out_dir, f"{r['accession']}.fna")
        if os.path.exists(out) and os.path.getsize(out) > 0:
            continue
        try:
            n = download_fasta(r["accession"], out)
            print(f"[genomes] {i}/{len(uniq)} {r['accession']} {r['organism']} ({n} file(s))")
        except Exception as e:  # noqa: BLE001
            failed += 1; print(f"[genomes] {r['accession']} FAILED ({e})", file=sys.stderr)
        time.sleep(0.34 if not os.environ.get("NCBI_API_KEY") else 0.1)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
