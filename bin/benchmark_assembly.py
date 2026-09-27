#!/usr/bin/env python3
"""Benchmark an assembly against a reference genome with minimap2 (W1.1).

Reports the QUAST-style numbers that matter for a fungal isolate assembly, without QUAST
(not in the images): aligned fraction of assembly and reference, NGA50 (aligned-block N50
against the reference length), alignment breakpoints per contig (misassembly candidates),
mismatches and indels per 100 kb from the cs tag, and unaligned contigs.

  benchmark_assembly.py --assembly asm.fa --reference ref.fa --out bench.json [--threads 8]
                        [--min-block 1000] [--paf existing.paf]

Breakpoint counting: a contig aligned in more than one block (on different reference
sequences, or on the same sequence but not collinear within --gap-tol) counts each extra
block as one breakpoint. Blocks shorter than --min-block are ignored for NGA50 and
breakpoints but still count toward the aligned fraction.
"""
from __future__ import annotations
import argparse, collections, itertools, json, os, re, shutil, subprocess, sys

CS_RE = re.compile(r"(:[0-9]+|\*[a-z][a-z]|[+-][a-z]+)")


def read_fasta_lengths(path):
    lens, name, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens


def run_minimap2(reference, assembly, threads):
    cmd = ["minimap2", "-cx", "asm5", "--cs", "-t", str(threads), reference, assembly]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"minimap2 failed: {out.stderr[-500:]}")
    return out.stdout.splitlines()


def parse_paf(lines):
    recs = []
    for line in lines:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        rec = {"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4],
               "t": f[5], "tlen": int(f[6]), "ts": int(f[7]), "te": int(f[8]), "matches": int(f[9]),
               "alen": int(f[10]), "mapq": int(f[11]), "tp": "P", "cs": None, "nm": None}
        for tag in f[12:]:
            if tag.startswith("tp:A:"):
                rec["tp"] = tag[5:]
            elif tag.startswith("cs:Z:"):
                rec["cs"] = tag[5:]
            elif tag.startswith("NM:i:"):
                rec["nm"] = int(tag[5:])
        recs.append(rec)
    return recs


def count_cs(cs):
    """(mismatches, insertions, deletions, inserted_bp, deleted_bp) from a cs string."""
    mm = ins = dele = ibp = dbp = 0
    if not cs:
        return mm, ins, dele, ibp, dbp
    for tok in CS_RE.findall(cs):
        if tok.startswith("*"):
            mm += 1
        elif tok.startswith("+"):
            ins += 1; ibp += len(tok) - 1
        elif tok.startswith("-"):
            dele += 1; dbp += len(tok) - 1
    return mm, ins, dele, ibp, dbp


def nga50(block_lengths, ref_total):
    acc = 0
    for L in sorted(block_lengths, reverse=True):
        acc += L
        if acc >= ref_total / 2:
            return L
    return 0


def merge_intervals(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return sum(e - s for s, e in out)


def benchmark(asm_lens, ref_lens, recs, min_block=1000, gap_tol=5000):
    primary = [r for r in recs if r["tp"] == "P"]
    asm_total, ref_total = sum(asm_lens.values()), sum(ref_lens.values())
    # coverage of assembly and reference by primary alignments (overlaps merged)
    q_iv, t_iv = collections.defaultdict(list), collections.defaultdict(list)
    for r in primary:
        q_iv[r["q"]].append((r["qs"], r["qe"])); t_iv[r["t"]].append((r["ts"], r["te"]))
    asm_aligned = sum(merge_intervals(v) for v in q_iv.values())
    ref_aligned = sum(merge_intervals(v) for v in t_iv.values())
    # per-contig breakpoints: blocks >= min_block, sorted by query start; a new block that is
    # on another reference sequence, another strand, or not collinear (gap > gap_tol) breaks
    breakpoints, blocks_for_nga = 0, []
    per_contig = {}
    for q, rs in itertools.groupby(sorted(primary, key=lambda r: (r["q"], r["qs"])), key=lambda r: r["q"]):
        rs = [r for r in rs if r["alen"] >= min_block]
        bp = 0
        for prev, cur in zip(rs, rs[1:]):
            collinear = (prev["t"] == cur["t"] and prev["strand"] == cur["strand"] and
                         (abs(cur["ts"] - prev["te"]) <= gap_tol if cur["strand"] == "+" else abs(prev["ts"] - cur["te"]) <= gap_tol))
            if not collinear:
                bp += 1
        breakpoints += bp
        blocks_for_nga += [r["te"] - r["ts"] for r in rs]
        per_contig[q] = {"length": asm_lens.get(q), "blocks": len(rs), "breakpoints": bp,
                         "targets": sorted({r["t"] for r in rs})}
    mm = ins = dele = ibp = dbp = 0
    for r in primary:
        a, b, c, d, e = count_cs(r["cs"])
        mm += a; ins += b; dele += c; ibp += d; dbp += e
    per100kb = (lambda n: round(100000 * n / ref_aligned, 2) if ref_aligned else None)
    unaligned = [q for q in asm_lens if q not in q_iv]
    return {
        "assembly_len": asm_total, "n_contigs": len(asm_lens), "reference_len": ref_total, "n_reference_seqs": len(ref_lens),
        "assembly_aligned_frac": round(asm_aligned / asm_total, 4) if asm_total else None,
        "reference_covered_frac": round(ref_aligned / ref_total, 4) if ref_total else None,
        "nga50": nga50(blocks_for_nga, ref_total),
        "n_alignment_blocks": len([r for r in primary if r["alen"] >= min_block]),
        "breakpoints": breakpoints,
        "contigs_with_breakpoints": sum(1 for c in per_contig.values() if c["breakpoints"]),
        "mismatches": mm, "indels": ins + dele, "inserted_bp": ibp, "deleted_bp": dbp,
        "mismatches_per_100kb": per100kb(mm), "indels_per_100kb": per100kb(ins + dele),
        "unaligned_contigs": len(unaligned), "unaligned_bp": sum(asm_lens[q] for q in unaligned),
        "per_contig": per_contig,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True); ap.add_argument("--reference", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-block", type=int, default=1000); ap.add_argument("--gap-tol", type=int, default=5000)
    ap.add_argument("--paf", default=None, help="use an existing minimap2 PAF (with --cs) instead of running minimap2")
    ap.add_argument("--sample", default=None)
    a = ap.parse_args()
    if a.paf:
        with open(a.paf) as fh:
            lines = fh.read().splitlines()
    else:
        if not shutil.which("minimap2"):
            sys.exit("[benchmark_assembly] minimap2 not found on PATH")
        lines = run_minimap2(a.reference, a.assembly, a.threads)
    res = benchmark(read_fasta_lengths(a.assembly), read_fasta_lengths(a.reference), parse_paf(lines), a.min_block, a.gap_tol)
    res = {"sample": a.sample or os.path.basename(a.assembly), "assembly": os.path.abspath(a.assembly),
           "reference": os.path.abspath(a.reference), **res}
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"[benchmark_assembly] {res['sample']}: {res['assembly_len']/1e6:.2f} Mb in {res['n_contigs']} contigs; "
          f"aligned {100*res['assembly_aligned_frac']:.1f} % of assembly / {100*res['reference_covered_frac']:.1f} % of reference; "
          f"NGA50 {res['nga50']:,}; breakpoints {res['breakpoints']}; mismatches/100kb {res['mismatches_per_100kb']}; "
          f"indels/100kb {res['indels_per_100kb']}; unaligned contigs {res['unaligned_contigs']}")


if __name__ == "__main__":
    main()
