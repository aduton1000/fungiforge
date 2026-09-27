#!/usr/bin/env python3
"""Cohort phylogenomics and clonality (W3.1, run-level stage 16).

  cohort_phylo.py --outdir cohort --threads N --assemblies a.fa b.fa ... --busco a.faa b.faa ...
                  --species a.species.txt ... --qc a.assemblyqc.json ...
                  [--min-frac 0.8] [--min-isolates 4] [--tree iqtree|fasttree] [--snp-threshold 100] [--ani-cluster 95]

1. skani triangle all-vs-all -> ANI matrix; species clusters by single linkage at --ani-cluster.
2. BUSCO single-copy orthologs present in >= --min-frac of the isolates -> MAFFT per gene ->
   gap-trimmed (columns with > 50 % gaps dropped) -> concatenated -> IQ-TREE 3 (LG+G, 1000
   ultrafast bootstraps) or FastTree; skipped below --min-isolates.
3. Within each cluster of >= 2 isolates: the most complete assembly (BUSCO, then N50) is the
   reference; every member is aligned to it (minimap2 asm5, --cs) and paftools.js call gives
   aligned regions (R) and variants (V); core = reference positions aligned in every member;
   pairwise SNP distance = SNVs differing between two members inside the core (indels ignored:
   assembly indels are the ONT/homopolymer error mode); clonal groups by single linkage at
   --snp-threshold. Outputs: ani_matrix.tsv, species_clusters.tsv, cohort.treefile,
   snp_distances.<cluster>.tsv, clonal_groups.tsv, cohort.json."""
from __future__ import annotations
import argparse, collections, json, os, shutil, subprocess, sys


def sample_of(path):
    base = os.path.basename(path)
    for suf in (".nuclear.fasta", ".busco_sc.faa", ".species.txt", ".assemblyqc.json", ".fasta", ".fa", ".faa", ".fna"):
        if base.endswith(suf):
            return base[:-len(suf)]
    return base.split(".")[0]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    if not path or not os.path.exists(path):
        return seqs
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def run(cmd, log=None, **kw):
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if log:
        with open(log, "a") as fh:
            fh.write(f"$ {' '.join(map(str, cmd))}\n{p.stderr}\n")
    return p


# ---------- ANI clusters ----------
def parse_skani_edges(text):
    """skani triangle -E output -> {(a, b): (ani, af_max)} keyed by sample names."""
    edges = {}
    lines = text.splitlines()
    if not lines:
        return edges
    header = lines[0].split("\t"); ci = {c: i for i, c in enumerate(header)}
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) < 3:
            continue
        try:
            ani = float(f[ci.get("ANI", 2)]); af = max(float(f[ci.get("Align_fraction_ref", 3)]), float(f[ci.get("Align_fraction_query", 4)]))
        except (ValueError, IndexError):
            continue
        a, b = sample_of(f[ci.get("Ref_file", 0)]), sample_of(f[ci.get("Query_file", 1)])
        edges[tuple(sorted((a, b)))] = (ani, af)
    return edges


def single_linkage(items, linked):
    """linked(a, b) -> bool; returns list of clusters (sorted lists)."""
    parent = {x: x for x in items}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    items = list(items)
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if linked(items[i], items[j]):
                parent[find(items[i])] = find(items[j])
    groups = collections.defaultdict(list)
    for x in items:
        groups[find(x)].append(x)
    return sorted((sorted(g) for g in groups.values()), key=lambda g: (-len(g), g[0]))


def ani_clusters(samples, edges, threshold=95.0, min_af=15.0):
    return single_linkage(samples, lambda a, b: edges.get(tuple(sorted((a, b))), (0, 0))[0] >= threshold and edges.get(tuple(sorted((a, b))), (0, 0))[1] >= min_af)


# ---------- SNP distances ----------
def parse_paftools_call(text):
    """R lines -> aligned intervals (0-based half-open) per ref contig; V lines -> SNVs {(contig, pos0): alt}."""
    regions, snvs, n_indel = collections.defaultdict(list), {}, 0
    for line in text.splitlines():
        f = line.split("\t")
        if f[0] == "R" and len(f) >= 4:
            regions[f[1]].append((int(f[2]), int(f[3])))
        elif f[0] == "V" and len(f) >= 8:
            ref, alt = f[6], f[7]
            if len(ref) == 1 and len(alt) == 1 and "-" not in (ref, alt):
                snvs[(f[1], int(f[2]))] = alt.upper()
            else:
                n_indel += 1
    return regions, snvs, n_indel


def intersect_regions(region_sets):
    """Intersection of per-member aligned intervals per contig -> {contig: [(s, e)]}."""
    if not region_sets:
        return {}
    contigs = set.intersection(*(set(r) for r in region_sets))
    core = {}
    for c in contigs:
        current = sorted(region_sets[0][c])
        for rs in region_sets[1:]:
            other = sorted(rs[c]); out = []; i = j = 0
            while i < len(current) and j < len(other):
                s, e = max(current[i][0], other[j][0]), min(current[i][1], other[j][1])
                if s < e:
                    out.append((s, e))
                if current[i][1] < other[j][1]:
                    i += 1
                else:
                    j += 1
            current = out
        core[c] = current
    return core


def in_core(core, contig, pos):
    for s, e in core.get(contig, []):
        if s <= pos < e:
            return True
    return False


def snp_distances(members, calls, ref_name):
    """members: sample names (incl. ref); calls: {sample: (regions, snvs)}; ref has no call (all-ref).
    Returns (core_bp, matrix {a: {b: dist}}) using SNVs inside the core."""
    region_sets = [calls[m][0] for m in members if m != ref_name]
    core = intersect_regions(region_sets) if region_sets else {}
    core_bp = sum(e - s for iv in core.values() for s, e in iv)
    alleles = {m: {k: v for k, v in calls[m][1].items() if in_core(core, *k)} for m in members if m != ref_name}
    alleles[ref_name] = {}
    dist = {a: {} for a in members}
    for i, a in enumerate(members):
        for b in members[i:]:
            keys = set(alleles[a]) | set(alleles[b])
            d = sum(1 for k in keys if alleles[a].get(k, "REF") != alleles[b].get(k, "REF"))
            dist[a][b] = dist[b][a] = d
    return core_bp, dist


def choose_reference(members, qc):
    """Most complete assembly: highest BUSCO, then N50; qc = {sample: assemblyqc dict}."""
    def key(m):
        q = qc.get(m, {}) or {}
        return (q.get("busco_complete") or 0, q.get("n50") or 0)
    return max(members, key=key)


# ---------- alignment trimming / concatenation ----------
def trim_columns(aln, max_gap_frac=0.5):
    """aln: {name: aligned seq}; drop columns where more than max_gap_frac of sequences have a gap."""
    if not aln:
        return aln
    names = list(aln); L = len(next(iter(aln.values()))); n = len(names)
    keep = [i for i in range(L) if sum(1 for nm in names if aln[nm][i] in "-?X") <= max_gap_frac * n]
    return {nm: "".join(aln[nm][i] for i in keep) for nm in names}


def concatenate(gene_alns, samples):
    """gene_alns: [{sample: seq}] -> supermatrix {sample: seq} (missing gene = gaps) and partitions."""
    parts, cat, pos = [], {s: [] for s in samples}, 1
    for gi, aln in enumerate(gene_alns):
        L = len(next(iter(aln.values())))
        for s in samples:
            cat[s].append(aln.get(s, "-" * L))
        parts.append((f"gene{gi + 1}", pos, pos + L - 1)); pos += L
    return {s: "".join(v) for s, v in cat.items()}, parts


def select_genes(busco_sets, samples, min_frac=0.8):
    """busco_sets: {sample: {busco_id: seq}} -> BUSCO ids single-copy in >= min_frac of samples."""
    counts = collections.Counter(b for s in samples for b in busco_sets.get(s, {}))
    need = max(2, int(round(min_frac * len(samples))))
    return sorted(b for b, c in counts.items() if c >= need)


# ---------- driver ----------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True); ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--assemblies", nargs="+", required=True); ap.add_argument("--busco", nargs="*", default=[])
    ap.add_argument("--species", nargs="*", default=[]); ap.add_argument("--qc", nargs="*", default=[])
    ap.add_argument("--min-frac", type=float, default=0.8); ap.add_argument("--min-isolates", type=int, default=4)
    ap.add_argument("--tree", choices=["iqtree", "fasttree"], default="iqtree"); ap.add_argument("--snp-threshold", type=int, default=100)
    ap.add_argument("--ani-cluster", type=float, default=95.0); ap.add_argument("--max-genes", type=int, default=500)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True); log = os.path.join(a.outdir, "cohort.log")
    asm = {sample_of(p): p for p in a.assemblies}
    samples = sorted(asm)
    species = {}
    for p in a.species:
        try:
            species[sample_of(p)] = open(p).readline().split("\t")[0].strip()
        except Exception:  # noqa: BLE001
            pass
    qc = {}
    for p in a.qc:
        try:
            qc[sample_of(p)] = json.load(open(p))
        except Exception:  # noqa: BLE001
            pass
    doc = {"stage": "cohort", "n_isolates": len(samples), "isolates": samples, "species": species, "tools": {}}
    # 1. ANI
    edges = {}
    if len(samples) >= 2 and shutil.which("skani"):
        lst = os.path.join(a.outdir, "assemblies.txt"); open(lst, "w").write("\n".join(asm[s] for s in samples) + "\n")
        p = run(["skani", "triangle", "-l", lst, "-E", "-t", str(a.threads), "--min-af", "0"], log)
        if p.returncode == 0:
            edges = parse_skani_edges(p.stdout); doc["tools"]["skani"] = "ok"
            with open(os.path.join(a.outdir, "ani_matrix.tsv"), "w") as fh:
                fh.write("sample\t" + "\t".join(samples) + "\n")
                for x in samples:
                    # pairs skani never reports (too distant to screen) are NA, not 0
                    fh.write(x + "\t" + "\t".join("100" if x == y else (f"{edges[tuple(sorted((x, y)))][0]:.2f}" if tuple(sorted((x, y))) in edges else "NA") for y in samples) + "\n")
        else:
            doc["tools"]["skani"] = "failed"
    clusters = ani_clusters(samples, edges, a.ani_cluster) if edges else [[s] for s in samples]
    doc["species_clusters"] = [{"cluster": f"C{i + 1}", "members": c, "species": sorted({species.get(m, "unknown") for m in c})} for i, c in enumerate(clusters)]
    with open(os.path.join(a.outdir, "species_clusters.tsv"), "w") as fh:
        fh.write("cluster\tsample\tspecies\n")
        for c in doc["species_clusters"]:
            for m in c["members"]:
                fh.write(f"{c['cluster']}\t{m}\t{species.get(m, 'unknown')}\n")
    # 2. phylogenomics
    busco_sets = {sample_of(p): read_fasta(p) for p in a.busco}
    busco_sets = {s: {n.split("|")[0]: q for n, q in v.items()} for s, v in busco_sets.items()}
    with_busco = [s for s in samples if busco_sets.get(s)]
    doc["phylogenomics"] = {"n_isolates_with_busco": len(with_busco), "n_genes": 0, "tree": None}
    if len(with_busco) >= a.min_isolates and shutil.which("mafft"):
        genes = select_genes(busco_sets, with_busco, a.min_frac)[:a.max_genes]
        gdir = os.path.join(a.outdir, "genes"); os.makedirs(gdir, exist_ok=True)
        alns = []
        for g in genes:
            fa = os.path.join(gdir, f"{g}.faa")
            with open(fa, "w") as fh:
                for s in with_busco:
                    if g in busco_sets[s]:
                        fh.write(f">{s}\n{busco_sets[s][g]}\n")
            p = run(["mafft", "--auto", "--quiet", "--thread", str(a.threads), fa], log)
            if p.returncode != 0:
                continue
            aln = {}; name = None
            for line in p.stdout.splitlines():
                if line.startswith(">"):
                    name = line[1:].split()[0]; aln[name] = ""
                elif name:
                    aln[name] += line.strip()
            alns.append(trim_columns(aln))
        if alns:
            supermatrix, parts = concatenate(alns, with_busco)
            sm = os.path.join(a.outdir, "cohort_alignment.faa")
            with open(sm, "w") as fh:
                for s in with_busco:
                    fh.write(f">{s}\n{supermatrix[s]}\n")
            with open(os.path.join(a.outdir, "cohort_partitions.txt"), "w") as fh:
                for name, s, e in parts:
                    fh.write(f"LG, {name} = {s}-{e}\n")
            doc["phylogenomics"].update(n_genes=len(alns), alignment_length=len(next(iter(supermatrix.values()))), genes=genes[:len(alns)])
            tree = os.path.join(a.outdir, "cohort.treefile")
            if a.tree == "iqtree" and shutil.which("iqtree3"):
                p = run(["iqtree3", "-s", sm, "-m", "LG+G", "-B", "1000", "-T", str(a.threads), "--prefix", os.path.join(a.outdir, "cohort"), "-quiet", "-redo"], log)
                doc["tools"]["iqtree3"] = "ok" if p.returncode == 0 else "failed"
            elif shutil.which("FastTree"):
                p = run(["FastTree", "-lg", "-gamma", sm], log)
                if p.returncode == 0:
                    open(tree, "w").write(p.stdout); doc["tools"]["FastTree"] = "ok"
            if os.path.exists(tree):
                doc["phylogenomics"]["tree"] = os.path.basename(tree); doc["phylogenomics"]["newick"] = open(tree).read().strip()[:20000]
    else:
        doc["phylogenomics"]["note"] = f"tree needs >= {a.min_isolates} isolates with single-copy BUSCO proteins (have {len(with_busco)})"
    # 3. SNP distances within clusters
    doc["clusters_snp"] = []; clonal_rows = []
    for ci, c in enumerate(clusters):
        if len(c) < 2 or not shutil.which("minimap2") or not shutil.which("paftools.js"):
            continue
        ref = choose_reference(c, qc); calls = {}
        for m in c:
            if m == ref:
                continue
            p = run(["minimap2", "-cx", "asm5", "--cs", "-t", str(a.threads), asm[ref], asm[m]], log)
            if p.returncode != 0:
                continue
            paf = os.path.join(a.outdir, f"{m}_vs_{ref}.paf"); open(paf, "w").write(p.stdout)
            srt = subprocess.run(["sort", "-k6,6", "-k8,8n", paf], capture_output=True, text=True)
            q = run(["paftools.js", "call", "-l", "10000", "-L", "10000", "-"], log, input=srt.stdout)
            if q.returncode == 0:
                regions, snvs, n_indel = parse_paftools_call(q.stdout); calls[m] = (regions, snvs)
        members = [ref] + [m for m in c if m != ref and m in calls]
        if len(members) < 2:
            continue
        core_bp, dist = snp_distances(members, calls, ref)
        name = f"C{ci + 1}"
        with open(os.path.join(a.outdir, f"snp_distances.{name}.tsv"), "w") as fh:
            fh.write("sample\t" + "\t".join(members) + "\n")
            for x in members:
                fh.write(x + "\t" + "\t".join(str(dist[x][y]) for y in members) + "\n")
        groups = single_linkage(members, lambda x, y: dist[x][y] <= a.snp_threshold)
        for gi, g in enumerate(groups):
            for m in g:
                clonal_rows.append((name, f"{name}.G{gi + 1}", m, len(g)))
        doc["clusters_snp"].append({"cluster": name, "reference": ref, "members": members, "core_bp": core_bp,
                                    "pairwise": {f"{x}|{y}": dist[x][y] for i, x in enumerate(members) for y in members[i + 1:]},
                                    "clonal_groups": [g for g in groups if len(g) > 1], "snp_threshold": a.snp_threshold,
                                    "note": "assembly-based SNVs inside the core aligned to the cluster reference; indels excluded; ONT-only assemblies carry residual error"})
    with open(os.path.join(a.outdir, "clonal_groups.tsv"), "w") as fh:
        fh.write("cluster\tclonal_group\tsample\tgroup_size\n")
        for r in clonal_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    json.dump(doc, open(os.path.join(a.outdir, "cohort.json"), "w"), indent=2)
    print(f"[cohort_phylo] {len(samples)} isolates, {len(clusters)} ANI clusters, tree={'yes' if doc['phylogenomics'].get('tree') else 'no'}, "
          f"SNP clusters={len(doc['clusters_snp'])}, clonal groups={sum(len(c['clonal_groups']) for c in doc['clusters_snp'])}")


if __name__ == "__main__":
    main()
