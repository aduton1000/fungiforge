#!/usr/bin/env python3
"""Build the antifungal-resistance REFERENCE SEQUENCES that bin/af_resistance.py aligns
against, from the FungAMR catalog. FungAMR ships as a mutation table (no sequences), but
its `ref_seq_uniprot_accession` column ties every mutation to a reference protein — so we
fetch those from UniProt and write:

  $DB/fungamr/reference_proteins.faa   headers ">GENE__ACC__Species" (gene-first, no spaces,
                                       so af_resistance's gene-keyword index picks them up)
  $DB/fungamr/fungamr_mutations.tsv    gene, accession, species, mutation, drug, confidence

Run on the HOST (needs network):  python3 bin/build_fungamr_refs.py --data-dir "$DB"
"""
from __future__ import annotations
import argparse, csv, json, os, re, sys, time, urllib.request, urllib.parse, urllib.error
from collections import Counter, defaultdict

UNIPROT = "https://rest.uniprot.org/uniprotkb/accessions"


def sanitize(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", (s or "").strip()).strip("_")


def _get(batch):
    """One UniProt request for a batch; returns fasta text, or None on HTTP error."""
    url = f"{UNIPROT}?accessions={urllib.parse.quote(','.join(batch))}&format=fasta"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code == 400:
                return None          # a bad/obsolete accession in the batch -> caller splits
            time.sleep(3)
        except Exception:
            time.sleep(3)
    return None


def _parse(txt, out):
    name, buf = None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if name: out[name] = "".join(buf)
            m = re.match(r">\w+\|([^|]+)\|", line)   # >sp|ACC|NAME or >tr|ACC|NAME
            name = m.group(1) if m else line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name: out[name] = "".join(buf)


def fetch_fasta(accessions, chunk=80):
    """Fetch sequences, recursively splitting any batch that 400s so a single bad/obsolete
    accession never sinks its neighbours."""
    out, bad = {}, []
    def do(batch):
        if not batch:
            return
        txt = _get(batch)
        if txt is not None:
            _parse(txt, out)
        elif len(batch) == 1:
            bad.append(batch[0])
        else:
            mid = len(batch) // 2
            do(batch[:mid]); do(batch[mid:])
    for i in range(0, len(accessions), chunk):
        do(accessions[i:i + chunk])
        sys.stderr.write(f"[fetch] {min(i+chunk,len(accessions))}/{len(accessions)} requested -> {len(out)} seqs, {len(bad)} unresolved\n")
    if bad:
        sys.stderr.write(f"[fetch] {len(bad)} accessions unresolved (obsolete/merged): {','.join(bad[:10])}...\n")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--min-confidence", type=int, default=8,
                    help="keep mutations with positive confidence score <= this (1=strongest..8=clinical)")
    a = ap.parse_args()
    fdir = os.path.join(a.data_dir, "fungamr")
    tsv = os.path.join(fdir, "FungAMR_070425.tsv")
    if not os.path.exists(tsv):
        sys.exit(f"FungAMR table not found: {tsv}")

    rows = list(csv.DictReader(open(tsv, encoding="utf-8", errors="replace"), delimiter="\t"))
    acc_gene = defaultdict(Counter)      # accession -> Counter(gene)
    acc_species = defaultdict(Counter)   # accession -> Counter(species)
    muts = []                            # (gene, acc, species, mutation, drug, confidence)
    accessions = set()
    for r in rows:
        raw = (r.get("ref_seq_uniprot_accession") or "").strip()
        gene = (r.get("gene or protein name") or "").strip()
        sp = (r.get("species") or "").strip()
        mut = (r.get("mutation") or "").strip()
        drug = (r.get("drug") or "").strip()
        conf = (r.get("confidence score") or "").strip()
        if not raw or raw.lower() in ("nan", "na", "none", ""):
            continue
        # keep positively-scored resistance mutations up to the requested confidence tier
        try:
            cv = float(conf)
        except ValueError:
            cv = None
        for acc in re.split(r"[,; ]+", raw):
            acc = acc.strip()
            if not acc:
                continue
            accessions.add(acc)
            if gene: acc_gene[acc][gene] += 1
            if sp: acc_species[acc][sp] += 1
            if mut and cv is not None and 0 < cv <= a.min_confidence:
                muts.append((gene, acc, sp, mut, drug, conf))

    accessions = sorted(accessions)
    sys.stderr.write(f"[build] {len(rows)} rows -> {len(accessions)} unique accessions, "
                     f"{len(muts)} confidence-filtered mutations\n")

    seqs = fetch_fasta(accessions)

    faa = os.path.join(fdir, "reference_proteins.faa")
    n = 0
    with open(faa, "w") as out:
        for acc, seq in seqs.items():
            if not seq:
                continue
            gene = acc_gene[acc].most_common(1)[0][0] if acc_gene[acc] else "unknown"
            sp = acc_species[acc].most_common(1)[0][0] if acc_species[acc] else "unknown"
            out.write(f">{sanitize(gene)}__{acc}__{sanitize(sp)}\n{seq}\n")
            n += 1

    mtsv = os.path.join(fdir, "fungamr_mutations.tsv")
    with open(mtsv, "w") as out:
        out.write("gene\taccession\tspecies\tmutation\tdrug\tconfidence\n")
        for row in muts:
            out.write("\t".join(row) + "\n")

    json.dump({"accessions_requested": len(accessions), "sequences_written": n,
               "mutations": len(muts)}, open(os.path.join(fdir, "reference_build.json"), "w"), indent=2)
    print(f"[build_fungamr_refs] {n} reference proteins -> {faa}")
    print(f"[build_fungamr_refs] {len(muts)} mutations -> {mtsv}")


if __name__ == "__main__":
    main()
