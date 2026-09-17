#!/usr/bin/env python3
"""Cohort gene-cluster families (W3.2, run-level stage 16).

  bgc_families.py --outdir cohort --threads N --regions a.regions.gbk b.regions.gbk ... [--bgc a.bgc.json ...]
                  [--min-sim 0.5] [--min-pident 50] [--min-regions 2]

Two methods, chosen by --method (default auto):

  bigscape        BiG-SCAPE 2 `cluster` over the per-region GenBanks with a Pfam database
                  (--pfam-path, by default the funannotate database's Pfam-A.hmm): profile-HMM
                  domains, per-class distance and affinity-propagation GCFs at --gcf-cutoff.
                  Used automatically when `bigscape` and a Pfam file are both available.
  shared_proteins DIAMOND all-vs-all of the region proteins; two regions join a family when at
                  least --min-sim of the smaller region's proteins have a homologue (>= --min-pident
                  identity, >= 50 % query coverage) in the other, single linkage. The fallback when
                  BiG-SCAPE or Pfam is unavailable.

Each family lists its members per isolate, product types and the best MIBiG description (from
bgc.json, when KnownClusterBlast ran)."""
from __future__ import annotations
import argparse, collections, json, os, re, shutil, subprocess


def sample_of(path):
    base = os.path.basename(path)
    for suf in (".regions.gbk", ".bgc.json", ".gbk"):
        if base.endswith(suf):
            return base[:-len(suf)]
    return base.split(".")[0]


def split_regions(gbk_path, sample, out_dir):
    """Write each LOCUS..// block of a concatenated regions GenBank to its own file named
    <sample>_<record>.regionNNN.gbk (BiG-SCAPE selects input files by the string 'region')."""
    os.makedirs(out_dir, exist_ok=True)
    written, buf, locus = [], [], None
    with open(gbk_path, errors="replace") as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                buf, locus = [line], line.split()[1]
            elif buf is not None:
                buf.append(line)
                if line.startswith("//"):
                    n = len(written) + 1
                    out = os.path.join(out_dir, f"{sample}_{locus}.region{n:03d}.gbk")
                    with open(out, "w") as o:
                        o.writelines(buf)
                    written.append(out); buf, locus = None, None
    return written


def parse_bigscape_clustering(bs_dir):
    """All *_clustering_c*.tsv of a BiG-SCAPE run -> {gbk_stem: family_id} (class-prefixed ids)."""
    import glob as _glob
    fam = {}
    for f in sorted(_glob.glob(os.path.join(bs_dir, "**", "*_clustering_c*.tsv"), recursive=True)):
        cls = os.path.basename(os.path.dirname(f))
        with open(f) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            ci = {c: i for i, c in enumerate(header)}
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) <= max(ci.get("GBK", 1), ci.get("Family", 5)):
                    continue
                fam[p[ci.get("GBK", 1)]] = f"{cls}:{p[ci.get('Family', 5)]}"
    return fam


def run_bigscape(region_files, out_dir, pfam_path, threads=4, cutoff=0.3, log=None):
    """Run BiG-SCAPE 2 cluster; returns {gbk_stem: family} or None when it could not run."""
    if not shutil.which("bigscape") or not pfam_path or not os.path.exists(pfam_path):
        return None
    gbk_dir = os.path.join(out_dir, "bigscape_input")
    os.makedirs(gbk_dir, exist_ok=True)
    for f in region_files:
        dest = os.path.join(gbk_dir, os.path.basename(f))
        if not os.path.exists(dest):
            shutil.copy(f, dest)
    # BiG-SCAPE presses the Pfam file in place; when the staged copy is not pressed and its
    # directory is read-only, press a local copy instead.
    pfam = pfam_path
    if not all(os.path.exists(pfam_path + ext) for ext in (".h3f", ".h3i", ".h3m", ".h3p")):
        if os.access(os.path.dirname(pfam_path) or ".", os.W_OK):
            subprocess.run(["hmmpress", "-f", pfam_path], capture_output=True, text=True)
        else:
            local = os.path.join(out_dir, "Pfam-A.hmm")
            if not os.path.exists(local):
                shutil.copy(pfam_path, local)
            subprocess.run(["hmmpress", "-f", local], capture_output=True, text=True)
            pfam = local
    bs_out = os.path.join(out_dir, "bigscape")
    shutil.rmtree(bs_out, ignore_errors=True)
    cmd = ["bigscape", "cluster", "-i", gbk_dir, "-p", pfam, "-o", bs_out, "-c", str(threads),
           "--no-trees", "--include-singletons", "--gcf-cutoffs", str(cutoff), "--quiet"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if log:
        with open(log, "a") as fh:
            fh.write(f"$ {' '.join(cmd)}\n{p.stdout[-3000:]}\n{p.stderr[-3000:]}\n")
    if p.returncode != 0:
        return None
    fam = parse_bigscape_clustering(bs_out)
    return fam or None


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
    ap.add_argument("--method", choices=["auto", "bigscape", "shared_proteins"], default="auto")
    ap.add_argument("--pfam-path", default=None, help="Pfam-A.hmm for BiG-SCAPE (default: <data-dir>/funannotate/Pfam-A.hmm)")
    ap.add_argument("--data-dir", default=""); ap.add_argument("--gcf-cutoff", type=float, default=0.3)
    ap.add_argument("--diamond-hits", default=None, help="pre-computed DIAMOND table (tests)")
    ap.add_argument("--bigscape-dir", default=None, help="pre-computed BiG-SCAPE output directory (tests)")
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
    doc = {"n_regions": len(regions), "n_isolates": len({r["sample"] for r in regions}), "method": None, "min_sim": a.min_sim}
    if len(regions) < a.min_regions or not any(r["proteins"] for r in regions):
        doc.update(families=[], note=f"fewer than {a.min_regions} regions with proteins")
        json.dump(doc, open(os.path.join(a.outdir, "bgc_families.json"), "w"), indent=2)
        open(os.path.join(a.outdir, "bgc_families.tsv"), "w").write("family\tregion\tsample\tproducts\tn_proteins\tmibig_best\n")
        print(f"[bgc_families] {len(regions)} regions: nothing to cluster"); return
    # --- BiG-SCAPE path (primary when available) ---
    pfam = a.pfam_path or (os.path.join(a.data_dir, "funannotate", "Pfam-A.hmm") if a.data_dir else None)
    bs_fam, split_files = None, []
    if a.method in ("auto", "bigscape") or a.bigscape_dir:
        for p in a.regions:
            if os.path.exists(p) and os.path.getsize(p) > 0:
                split_files += [(f, sample_of(p)) for f in split_regions(p, sample_of(p), os.path.join(a.outdir, "regions_split"))]
    if a.bigscape_dir:
        bs_fam = parse_bigscape_clustering(a.bigscape_dir) or None
    elif a.method in ("auto", "bigscape") and split_files:
        bs_fam = run_bigscape([f for f, _ in split_files], a.outdir, pfam, a.threads, a.gcf_cutoff, os.path.join(a.outdir, "cohort.log"))
    if bs_fam is None and a.method == "bigscape":
        doc.update(families=[], method="bigscape", note="BiG-SCAPE could not run (bigscape or the Pfam database is missing)")
        json.dump(doc, open(os.path.join(a.outdir, "bgc_families.json"), "w"), indent=2)
        open(os.path.join(a.outdir, "bgc_families.tsv"), "w").write("family\tregion\tsample\tproducts\tn_proteins\tmibig_best\n")
        print("[bgc_families] BiG-SCAPE unavailable"); return
    if bs_fam:
        # each split file holds exactly one region, so BiG-SCAPE's GBK column maps to it directly
        split_regs = []
        for f, sample in sorted(split_files):
            for r in parse_regions(f, sample):
                r["file_stem"] = os.path.basename(f)[:-4]
                split_regs.append(r)
        groups = collections.defaultdict(list)
        for idx, r in enumerate(split_regs):
            groups[bs_fam.get(r["file_stem"]) or f"__unassigned_{idx}"].append(idx)
        fams = sorted((sorted(g) for g in groups.values()), key=lambda g: (-len(g), g[0]))
        doc.update(method="bigscape", gcf_cutoff=a.gcf_cutoff, n_regions=len(split_regs),
                   n_isolates=len({r["sample"] for r in split_regs}),
                   n_regions_assigned=sum(1 for r in split_regs if r["file_stem"] in bs_fam))
        write_families(a.outdir, split_regs, fams, mibig, doc)
        return
    doc["method"] = "diamond_shared_proteins"
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
    write_families(a.outdir, regions, fams, mibig, doc)


def write_families(outdir, regions, fams, mibig, doc):
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
    with open(os.path.join(outdir, "bgc_families.tsv"), "w") as fh:
        fh.write("family\tregion\tsample\tproducts\tn_proteins\tmibig_best\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    shared = [f for f in fam_docs if f["n_isolates"] > 1]
    doc.update(n_families=len(fams), n_families_shared=len(shared), n_singletons=sum(1 for f in fam_docs if f["n_regions"] == 1), families=fam_docs,
               known_families=[f["family"] for f in fam_docs if f["mibig_best"]])
    json.dump(doc, open(os.path.join(outdir, "bgc_families.json"), "w"), indent=2)
    print(f"[bgc_families] {len(regions)} regions from {doc['n_isolates']} isolates -> {len(fams)} families "
          f"({len(shared)} shared by >1 isolate) by {doc.get('method')}")


if __name__ == "__main__":
    main()
