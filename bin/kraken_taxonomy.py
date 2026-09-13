#!/usr/bin/env python3
"""Kraken2 report utilities shared by Stage 04 (contig decontamination) and the
read-level triage (bin/triage_reads.sh).

A Kraken2 *report* encodes the whole taxonomy it touched (indentation = depth, one
rank code per line), so lineages can be reconstructed from the report alone — no
nodes.dmp / names.dmp needed. Every taxid that appears in the per-sequence output
(k2.out) also appears in the report, which is all we need to assign a domain.

Sub-commands
  summarize        k2.report [--json]        -> domain/kingdom composition + top species
  classify-contigs k2.report k2.out assembly.fasta --out-fasta kept.fasta --out-json d.json
                   [--keep-all]              -> drop Bacteria/Archaea/Viruses/human contigs,
                                               write composition (by bp) + verdict
Verdict rules (fraction of assigned bp for contigs; of reads for triage):
  non_fungal : bacteria+archaea+viruses >= 50 %   (or human >= 50 %)
  fungal     : fungi >= 10 %  and  bacteria+archaea+viruses < 20 %
  likely_fungal : unclassified >= 50 % and bacteria+archaea+viruses < 10 %
                  (the PlusPF-8 database holds few fungal genomes, so a genuine fungal
                  isolate often reads as mostly unclassified with a fungal minority)
  mixed      : anything else
"""
from __future__ import annotations
import argparse, json, sys

FUNGI = 4751; HUMAN = 9606; ROOT = 1; CELLULAR = 131567
DOMAINS = {2: "Bacteria", 2157: "Archaea", 2759: "Eukaryota", 10239: "Viruses"}


def parse_report(path):
    """taxid -> dict(name, rank, parent, depth, clade_pct, clade_n, taxon_n)."""
    nodes, stack = {}, []   # stack of (depth, taxid)
    for line in open(path):
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        pct, clade_n, taxon_n, rank, taxid, raw = f[0], f[1], f[2], f[3], f[4], f[5]
        name = raw.lstrip(" ")
        depth = (len(raw) - len(name)) // 2
        taxid = int(taxid)
        while stack and stack[-1][0] >= depth:
            stack.pop()
        parent = stack[-1][1] if stack else None
        nodes[taxid] = {"name": name, "rank": rank, "parent": parent, "depth": depth,
                        "clade_pct": float(pct), "clade_n": int(clade_n), "taxon_n": int(taxon_n)}
        stack.append((depth, taxid))
    return nodes


def lineage(nodes, taxid):
    out, t, guard = [], taxid, 0
    while t is not None and t in nodes and guard < 100:
        out.append(t); t = nodes[t]["parent"]; guard += 1
    return out


def group_of(nodes, taxid):
    """'unclassified' | 'human' | 'fungi' | 'other_eukaryota' | 'bacteria' | 'archaea' | 'viruses' | 'unresolved'."""
    if taxid == 0 or taxid not in nodes:
        return "unclassified"
    lin = lineage(nodes, taxid)
    if HUMAN in lin: return "human"
    if FUNGI in lin: return "fungi"
    for t in lin:
        if t in DOMAINS:
            d = DOMAINS[t]
            return {"Bacteria": "bacteria", "Archaea": "archaea", "Viruses": "viruses"}.get(d, "other_eukaryota")
    return "unresolved"   # root / cellular organisms only


def verdict(frac):
    """frac: dict group -> fraction (0..1) of reads or bp, incl. 'unclassified'."""
    prok = frac.get("bacteria", 0) + frac.get("archaea", 0) + frac.get("viruses", 0)
    if frac.get("human", 0) >= 0.5:
        return "human"
    if prok >= 0.5:
        return "non_fungal"
    if frac.get("fungi", 0) >= 0.10 and prok < 0.20:
        return "fungal"
    if frac.get("unclassified", 0) >= 0.5 and prok < 0.10:
        return "likely_fungal"
    return "mixed"


def top_species(nodes, k=3):
    sp = [(t, n) for t, n in nodes.items() if n["rank"] == "S"]
    sp.sort(key=lambda x: -x[1]["clade_pct"])
    return [{"taxid": t, "name": n["name"], "pct": n["clade_pct"], "group": group_of(nodes, t)} for t, n in sp[:k]]


def cmd_summarize(a):
    nodes = parse_report(a.report)
    total = sum(n["taxon_n"] for n in nodes.values())
    groups = {}
    for t, n in nodes.items():
        g = group_of(nodes, t)
        groups[g] = groups.get(g, 0) + n["taxon_n"]
    frac = {g: (c / total if total else 0.0) for g, c in groups.items()}
    res = {"total": total, "pct": {g: round(100 * v, 2) for g, v in frac.items()},
           "top_species": top_species(nodes), "verdict": verdict(frac)}
    if a.json:
        print(json.dumps(res, indent=2))
    else:
        p = res["pct"]; ts = res["top_species"][0] if res["top_species"] else {"name": "-", "pct": 0}
        print("\t".join(str(x) for x in [total, p.get("unclassified", 0), p.get("bacteria", 0), p.get("archaea", 0),
                                         p.get("viruses", 0), p.get("fungi", 0), p.get("other_eukaryota", 0),
                                         p.get("human", 0), ts["name"], ts["pct"], res["verdict"]]))


def read_fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None: yield name, "".join(buf)


def cmd_classify(a):
    nodes = parse_report(a.report)
    assign = {}
    for line in open(a.k2out):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            assign[f[1]] = int(f[2]) if f[0] == "C" else 0
    bp, kept_n, kept_bp, total_bp = {}, 0, 0, 0
    drop_groups = {"bacteria", "archaea", "viruses", "human"}
    contigs = list(read_fasta(a.assembly))
    for name, seq in contigs:
        g = group_of(nodes, assign.get(name, 0)); L = len(seq)
        bp[g] = bp.get(g, 0) + L; total_bp += L
    frac = {g: (v / total_bp if total_bp else 0.0) for g, v in bp.items()}
    v = verdict(frac)
    # A non-fungal isolate is not a contaminated fungus: filtering would leave next to
    # nothing. Keep the assembly whole and let the verdict/QC gate carry the message.
    keep_all = a.keep_all or v in ("non_fungal", "human")
    with open(a.out_fasta, "w") as fh:
        for name, seq in contigs:
            g = group_of(nodes, assign.get(name, 0))
            if keep_all or g not in drop_groups:
                fh.write(f">{name}\n"); [fh.write(seq[i:i + 80] + "\n") for i in range(0, len(seq), 80)]
                kept_n += 1; kept_bp += len(seq)
    res = {"sample": a.sample, "stage": "decontam", "nuclear": a.out_fasta, "mito": a.mito,
           "n_contigs_in": len(contigs), "bp_in": total_bp, "n_contigs_kept": kept_n, "bp_kept": kept_bp,
           "dropped_bp_pct": round(100 * (total_bp - kept_bp) / total_bp, 2) if total_bp else 0.0,
           "composition_bp_pct": {g: round(100 * f, 2) for g, f in frac.items()},
           "top_species": top_species(nodes), "verdict": v, "filtered": not keep_all,
           "note": ("non-fungal or human-dominated assembly kept whole (nothing to decontaminate); "
                    "downstream fungal stages are not meaningful for this isolate" if keep_all and not a.keep_all else
                    "bacterial/archaeal/viral/human contigs removed")}
    json.dump(res, open(a.out_json, "w"), indent=2)
    print(f"[decontam] {a.sample}: {v}; kept {kept_n}/{len(contigs)} contigs, {kept_bp}/{total_bp} bp; "
          f"composition {res['composition_bp_pct']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("summarize"); s.add_argument("report"); s.add_argument("--json", action="store_true"); s.set_defaults(f=cmd_summarize)
    c = sub.add_parser("classify-contigs"); c.add_argument("report"); c.add_argument("k2out"); c.add_argument("assembly")
    c.add_argument("--sample", default="sample"); c.add_argument("--out-fasta", required=True); c.add_argument("--out-json", required=True)
    c.add_argument("--mito", default=""); c.add_argument("--keep-all", action="store_true"); c.set_defaults(f=cmd_classify)
    a = ap.parse_args(); a.f(a)


if __name__ == "__main__":
    main()
