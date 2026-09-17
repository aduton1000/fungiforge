#!/usr/bin/env python3
"""Cohort gene-cluster families (W3.2, run-level stage 16).

  bgc_families.py --outdir cohort --threads N --regions a.regions.gbk b.regions.gbk ... [--bgc a.bgc.json ...]
                  [--min-sim 0.5] [--min-pident 50] [--min-regions 2]

Every antiSMASH region GenBank (all isolates) contributes its CDS proteins; a DIAMOND all-vs-all
of those proteins gives, for each pair of regions, the fraction of the smaller region's proteins
with a homologue (identity >= --min-pident, E <= 1e-10, alignment >= 50 % of the query) in the other
region; regions are linked into a gene-cluster family (GCF) when that fraction is >= --min-sim
(single linkage). Each family lists its members per isolate, product types and the best MIBiG
description (from bgc.json, when KnownClusterBlast ran). This is the in-image approximation of
BiG-SCAPE's domain-based clustering; when `bigscape` is on PATH and MIBiG is staged the stage
can run it instead (--bigscape)."""
from __future__ import annotations
import argparse, collections, json, os, re, shutil, subprocess


def sample_of(path):
    base = os.path.basename(path)
    for suf in (".regions.gbk", ".bgc.json", ".gbk"):
        if base.endswith(suf):
            return base[:-len(suf)]
    return base.split(".")[0]


def parse_regions(gbk_path, sample):
    """Yield {id, sample, record, region_number, product, proteins: {name: seq}, length} per region record."""
    regions, rec = [], None
    with open(gbk_path, errors="replace") as fh:
        lines = fh.read().split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("LOCUS"):
            rec = {"record": line.split()[1], "region_number": None, "product": [], "proteins": {}, "length": int(line.split()[2]) if line.split()[2].isdigit() else None}
        elif rec is not None and line.startswith("     region ") or (rec is not None and line.startswith("     region")):
            j = i + 1
            while j < len(lines) and lines[j].startswith("                     "):
                q = lines[j].strip()
                m = re.match(r'/region_number="?(\d+)"?', q)
                if m:
                    rec["region_number"] = int(m.group(1))
                m = re.match(r'/product="([^"]+)"', q)
                if m:
                    rec["product"].append(m.group(1))
                j += 1
        elif rec is not None and line.startswith("     CDS "):
            j = i + 1; name = None; tr = None
            while j < len(lines) and lines[j].startswith("                     "):
                q = lines[j].strip()
                m = re.match(r'/(locus_tag|protein_id|gene)="([^"]+)"', q)
                if m and name is None:
                    name = m.group(2)
                if q.startswith('/translation="'):
                    buf = q[len('/translation="'):]
                    while not buf.endswith('"'):
                        j += 1; buf += lines[j].strip()
                    tr = buf.rstrip('"')
                j += 1
            if tr:
                rec["proteins"][name or f"cds{len(rec['proteins']) + 1}"] = tr
        elif line.startswith("//") and rec is not None:
            rec["sample"] = sample
            rec["id"] = f"{sample}|{rec['record']}|r{rec['region_number'] if rec['region_number'] is not None else len(regions) + 1}"
            regions.append(rec); rec = None
        i += 1
    return regions


def parse_diamond(text, min_pident=50.0, min_cov=0.5):
    """outfmt 6 qseqid sseqid pident length qlen -> {query_protein: {subject_protein}} (self excluded)."""
    hits = collections.defaultdict(set)
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        try:
            pid, alen, qlen = float(f[2]), int(f[3]), int(f[4])
        except ValueError:
            continue
        if f[0] != f[1] and pid >= min_pident and alen >= min_cov * qlen:
            hits[f[0]].add(f[1])
    return hits


def region_similarity(regions, hits):
    """{(i, j): fraction of the smaller region's proteins with a homologue in the other}."""
    owner = {}
    for idx, r in enumerate(regions):
        for p in r["proteins"]:
            owner[f"{r['id']}||{p}"] = idx
    shared = collections.defaultdict(set)
    for q, subs in hits.items():
        qi = owner.get(q)
        if qi is None:
            continue
        for s in subs:
            si = owner.get(s)
            if si is not None and si != qi:
                shared[(qi, si)].add(q)
    sim = {}
    for (qi, si), qs in shared.items():
        key = tuple(sorted((qi, si)))
        n_small = min(len(regions[qi]["proteins"]), len(regions[si]["proteins"]))
        if n_small:
            frac = len(qs) / n_small
            sim[key] = max(sim.get(key, 0.0), frac)
    return sim


def families(regions, sim, min_sim=0.5):
    parent = list(range(len(regions)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (i, j), s in sim.items():
        if s >= min_sim:
            parent[find(i)] = find(j)
    groups = collections.defaultdict(list)
    for i in range(len(regions)):
        groups[find(i)].append(i)
    return sorted((sorted(g) for g in groups.values()), key=lambda g: (-len(g), g[0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True); ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--regions", nargs="*", default=[]); ap.add_argument("--bgc", nargs="*", default=[])
    ap.add_argument("--min-sim", type=float, default=0.5); ap.add_argument("--min-pident", type=float, default=50.0); ap.add_argument("--min-regions", type=int, default=2)
    ap.add_argument("--diamond-hits", default=None, help="pre-computed DIAMOND table (tests)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    regions = []
    for p in a.regions:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            regions.extend(parse_regions(p, sample_of(p)))
    mibig = {}
    for p in a.bgc:
        try:
            d = json.load(open(p)); s = sample_of(p)
            for key, best in (d.get("knownclusterblast", {}).get("best_hits") or {}).items():
                rec, rn = key.split("|")
                mibig[f"{s}|{rec}|r{rn}"] = best
        except Exception:  # noqa: BLE001
            pass
    doc = {"n_regions": len(regions), "n_isolates": len({r["sample"] for r in regions}), "method": "diamond_shared_proteins", "min_sim": a.min_sim}
    if len(regions) < a.min_regions or not any(r["proteins"] for r in regions):
        doc.update(families=[], note=f"fewer than {a.min_regions} regions with proteins")
        json.dump(doc, open(os.path.join(a.outdir, "bgc_families.json"), "w"), indent=2)
        open(os.path.join(a.outdir, "bgc_families.tsv"), "w").write("family\tregion\tsample\tproducts\tn_proteins\tmibig_best\n")
        print(f"[bgc_families] {len(regions)} regions: nothing to cluster"); return
    faa = os.path.join(a.outdir, "bgc_proteins.faa")
    with open(faa, "w") as fh:
        for r in regions:
            for name, seq in r["proteins"].items():
                fh.write(f">{r['id']}||{name}\n{seq}\n")
    text = None
    if a.diamond_hits and os.path.exists(a.diamond_hits):
        text = open(a.diamond_hits).read()
    elif shutil.which("diamond"):
        db = os.path.join(a.outdir, "bgc_proteins")
        subprocess.run(["diamond", "makedb", "--in", faa, "-d", db, "--quiet"], check=False)
        p = subprocess.run(["diamond", "blastp", "-q", faa, "-d", db, "-p", str(a.threads), "-e", "1e-10", "-k", "200", "--quiet", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen"],
                           capture_output=True, text=True)
        text = p.stdout if p.returncode == 0 else None
    if text is None:
        doc.update(families=[], note="DIAMOND not available: families not computed")
        json.dump(doc, open(os.path.join(a.outdir, "bgc_families.json"), "w"), indent=2)
        open(os.path.join(a.outdir, "bgc_families.tsv"), "w").write("family\tregion\tsample\tproducts\tn_proteins\tmibig_best\n")
        print("[bgc_families] DIAMOND unavailable"); return
    sim = region_similarity(regions, parse_diamond(text, a.min_pident))
    fams = families(regions, sim, a.min_sim)
    rows, fam_docs = [], []
    for fi, g in enumerate(fams):
        fid = f"GCF{fi + 1:04d}"
        members = [regions[i] for i in g]
        best = [mibig.get(m["id"]) for m in members]
        best = [b for b in best if b]
        top = max(best, key=lambda b: (b.get("core_gene_hits") or 0, b.get("blast_score") or 0)) if best else None
        fam_docs.append({"family": fid, "n_regions": len(members), "n_isolates": len({m["sample"] for m in members}),
                         "isolates": sorted({m["sample"] for m in members}), "products": sorted({p for m in members for p in m["product"]}),
                         "mibig_best": top, "members": [m["id"] for m in members]})
        for m in members:
            rows.append((fid, m["id"], m["sample"], ";".join(m["product"]) or "unknown", len(m["proteins"]), (mibig.get(m["id"]) or {}).get("description") or ""))
    with open(os.path.join(a.outdir, "bgc_families.tsv"), "w") as fh:
        fh.write("family\tregion\tsample\tproducts\tn_proteins\tmibig_best\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    shared = [f for f in fam_docs if f["n_isolates"] > 1]
    doc.update(n_families=len(fams), n_families_shared=len(shared), n_singletons=sum(1 for f in fam_docs if f["n_regions"] == 1), families=fam_docs,
               known_families=[f["family"] for f in fam_docs if f["mibig_best"]])
    json.dump(doc, open(os.path.join(a.outdir, "bgc_families.json"), "w"), indent=2)
    print(f"[bgc_families] {len(regions)} regions from {doc['n_isolates']} isolates -> {len(fams)} families ({len(shared)} shared by >1 isolate)")


if __name__ == "__main__":
    main()
