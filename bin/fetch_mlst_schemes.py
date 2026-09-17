#!/usr/bin/env python3
"""Stage the fungal PubMLST schemes for `mlst` (W2.3).

`mlst` ships bacterial schemes only; the fungal schemes are downloaded from the PubMLST REST
API into the layout mlst expects (--datadir): <out_dir>/pubmlst/<scheme>/<scheme>.txt (profiles)
and <scheme>/<locus>.tfa (alleles). The BLAST database mlst needs (--blastdb) is built inside
the pipeline task from these files (makeblastdb runs in the image, so no BLAST is needed here).

  fetch_mlst_schemes.py --out-dir <data_dir>/mlst [--schemes afumigatus,calbicans,cglabrata,ctropicalis,ckrusei]

Scheme names are the PubMLST database keys (pubmlst_<scheme>_seqdef), scheme id 1 in each.
"""
from __future__ import annotations
import argparse, json, os, sys, time, urllib.request

REST = "https://rest.pubmlst.org/db/pubmlst_{scheme}_seqdef/schemes/1"
DEFAULT = "afumigatus,calbicans,cglabrata,ctropicalis,ckrusei"


def _get(url, retries=5):
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "fungiforge"}), timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** i)
    raise RuntimeError(f"request failed: {url} ({last})")


def fetch_scheme(scheme, out_dir, log=print):
    meta = json.loads(_get(REST.format(scheme=scheme)))
    d = os.path.join(out_dir, "pubmlst", scheme)
    os.makedirs(d, exist_ok=True)
    loci = []
    for url in meta.get("loci", []):
        locus = url.rstrip("/").split("/")[-1]
        fa = _get(url.rstrip("/") + "/alleles_fasta")
        with open(os.path.join(d, f"{locus}.tfa"), "w") as fh:
            fh.write(fa if fa.endswith("\n") else fa + "\n")
        loci.append(locus)
    prof = _get(meta["profiles_csv"]) if meta.get("profiles_csv") else None
    if prof:
        with open(os.path.join(d, f"{scheme}.txt"), "w") as fh:
            fh.write(prof if prof.endswith("\n") else prof + "\n")
    n_prof = (prof.count("\n") - 1) if prof else 0
    log(f"[mlst] {scheme}: {len(loci)} loci ({', '.join(loci)}), {n_prof} profiles")
    return {"loci": loci, "profiles": n_prof, "description": meta.get("description")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--schemes", default=DEFAULT)
    a = ap.parse_args()
    summary, failed = {}, []
    for s in [x for x in a.schemes.split(",") if x]:
        try:
            summary[s] = fetch_scheme(s, a.out_dir)
        except Exception as e:  # noqa: BLE001
            print(f"[mlst] {s}: FAILED ({e})", file=sys.stderr); failed.append(s)
    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "mlst_manifest.json"), "w") as fh:
        json.dump({"date": time.strftime("%Y-%m-%d"), "schemes": summary, "failed": failed}, fh, indent=2)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
