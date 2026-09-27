#!/usr/bin/env python3
"""Verify a control isolate's genotype FROM ITS READS at a reference locus (W1.1).

A control is only a control if its reads carry the genotype claimed for it — the ENA reads
deposited under an isolate name can differ from the isolate whose assembly was published
(found with C87, 2026-09-15). This tool maps reads to a small wild-type locus and answers,
per check, whether the reads carry (a) an insertion of a given length range inside a window
(promoter tandem repeats such as TR34/TR46) or (b) an alternative base at a position (a
point mutation such as L98H). No assembly is involved, so a duplication an assembler might
collapse is still seen.

  control_genotype.py --spec test/controls/afum_cyp51A_controls.json --control TR34_L98H \
                      --reads reads.fq.gz [--reads more.fq.gz] [--platform ont|illumina] \
                      [--threads 8] [--out genotype.json] [--min-reads 10] [--min-frac 0.8]

The spec (JSON) names the locus FASTA (relative to the spec file), the checks, and the
named controls with the expected outcome of each check. --control picks which expectations
to apply; without it the observations are reported and no verdict is made.

Needs minimap2 on PATH. SAM parsing is pure Python (no pysam).
"""
from __future__ import annotations
import argparse, collections, gzip, json, os, re, shutil, subprocess, sys

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_cigar(cigar):
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def walk(pos, cigar_ops, seq):
    """Yield (ref_pos_1based, op, length, query_base_or_inserted_seq) for every CIGAR element.
    ref_pos is the reference position the element starts at (for I: the position it follows)."""
    r, q = pos, 0
    for n, op in cigar_ops:
        if op in ("M", "=", "X"):
            for i in range(n):
                yield r + i, "M", 1, seq[q + i] if q + i < len(seq) else "N"
            r += n; q += n
        elif op == "I":
            yield r - 1, "I", n, seq[q:q + n]
            q += n
        elif op == "D":
            yield r, "D", n, ""
            r += n
        elif op == "N":
            r += n
        elif op == "S":
            q += n
        # H, P consume nothing we track


def evaluate_reads(sam_iter, checks, min_mapq=1):
    """Consume SAM lines; return per-check observations."""
    obs = {}
    for c in checks:
        if c["type"] == "insertion":
            obs[c["name"]] = {"type": "insertion", "window": c["window"], "spanning": 0, "with_insertion": 0, "insertion_lengths": collections.Counter()}
        elif c["type"] == "snv":
            obs[c["name"]] = {"type": "snv", "pos": c["pos"], "ref": c["ref"].upper(), "alt": c["alt"].upper(), "depth": 0, "bases": collections.Counter()}
    for line in sam_iter:
        if not line or line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        flag, pos, mapq, cigar, seq = int(f[1]), int(f[3]), int(f[4]), f[5], f[9]
        if flag & 0x904 or cigar == "*" or mapq < min_mapq:      # unmapped, secondary, supplementary
            continue
        ops = parse_cigar(cigar)
        ref_len = sum(n for n, op in ops if op in "MDN=X")
        start, end = pos, pos + ref_len - 1
        ins_in_window = {c["name"]: [] for c in checks if c["type"] == "insertion"}
        bases_at = {}
        for rpos, op, n, payload in walk(pos, ops, seq):
            if op == "I":
                for c in checks:
                    if c["type"] == "insertion" and c["window"][0] - 5 <= rpos <= c["window"][1] + 5:
                        ins_in_window[c["name"]].append(n)
            elif op == "M":
                bases_at[rpos] = payload.upper()
            elif op == "D":
                for i in range(n):
                    bases_at[rpos + i] = "-"
        for c in checks:
            o = obs[c["name"]]
            if c["type"] == "insertion":
                w0, w1 = c["window"]
                if start <= w0 - 20 and end >= w1 + 20:          # read spans the window with margin
                    o["spanning"] += 1
                    lens = [L for L in ins_in_window[c["name"]] if c.get("min_len", 1) <= L <= c.get("max_len", 10**9)]
                    if lens:
                        o["with_insertion"] += 1
                        o["insertion_lengths"][max(lens)] += 1
            else:
                b = bases_at.get(c["pos"])
                if b is not None and start <= c["pos"] <= end:
                    o["depth"] += 1
                    o["bases"][b] += 1
    for o in obs.values():
        if o["type"] == "insertion":
            o["insertion_lengths"] = dict(sorted(o["insertion_lengths"].items()))
            o["frac"] = round(o["with_insertion"] / o["spanning"], 3) if o["spanning"] else None
        else:
            o["bases"] = dict(o["bases"])
            o["alt_frac"] = round(o["bases"].get(o["alt"], 0) / o["depth"], 3) if o["depth"] else None
            o["ref_frac"] = round(o["bases"].get(o["ref"], 0) / o["depth"], 3) if o["depth"] else None
    return obs


def call(obs, min_reads, min_frac):
    """Per check: present | absent | insufficient."""
    calls = {}
    for name, o in obs.items():
        n = o["spanning"] if o["type"] == "insertion" else o["depth"]
        frac = o["frac"] if o["type"] == "insertion" else o["alt_frac"]
        if n < min_reads or frac is None:
            calls[name] = "insufficient"
        elif frac >= min_frac:
            calls[name] = "present"
        elif frac <= 1 - min_frac:
            calls[name] = "absent"
        else:
            calls[name] = "mixed"
    return calls


def run_minimap2(locus, reads, platform, threads):
    preset = "map-ont" if platform == "ont" else "sr"
    cmd = ["minimap2", "-a", "-x", preset, "-t", str(threads), "--secondary=no", locus] + reads
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    for line in p.stdout:
        yield line
    p.wait()
    if p.returncode != 0:
        raise RuntimeError(f"minimap2 exited {p.returncode}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--control", default=None, help="named control in the spec whose expectations to apply")
    ap.add_argument("--reads", action="append", required=True)
    ap.add_argument("--platform", choices=["ont", "illumina"], default="ont")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--min-frac", type=float, default=0.8)
    ap.add_argument("--sam", default=None, help="parse this SAM instead of running minimap2 (tests)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    with open(a.spec) as fh:
        spec = json.load(fh)
    locus = os.path.join(os.path.dirname(os.path.abspath(a.spec)), spec["locus"])
    if a.sam:
        opener = gzip.open if a.sam.endswith(".gz") else open
        with opener(a.sam, "rt") as fh:
            obs = evaluate_reads(fh, spec["checks"])
    else:
        if not shutil.which("minimap2"):
            sys.exit("[control_genotype] minimap2 not found on PATH")
        obs = evaluate_reads(run_minimap2(locus, a.reads, a.platform, a.threads), spec["checks"])
    calls = call(obs, a.min_reads, a.min_frac)
    result = {"spec": os.path.abspath(a.spec), "locus": locus, "reads": a.reads, "platform": a.platform,
              "observations": obs, "calls": calls, "control": a.control}
    if a.control:
        exp = spec["controls"][a.control]["expect"]
        verdict = {}
        for name, want in exp.items():
            got = calls.get(name, "insufficient")
            verdict[name] = "PASS" if got == ("present" if want else "absent") else ("INSUFFICIENT" if got == "insufficient" else "FAIL")
        result["expected"] = exp
        result["verdict"] = verdict
        result["passed"] = all(v == "PASS" for v in verdict.values())
    for name, o in obs.items():
        if o["type"] == "insertion":
            print(f"[control_genotype] {name:8s} insertion {o['window']}: {o['with_insertion']}/{o['spanning']} spanning reads (frac {o['frac']}), lengths {o['insertion_lengths']} -> {calls[name]}")
        else:
            print(f"[control_genotype] {name:8s} snv pos {o['pos']} {o['ref']}>{o['alt']}: depth {o['depth']}, bases {o['bases']} (alt frac {o['alt_frac']}) -> {calls[name]}")
    if a.control:
        print(f"[control_genotype] control {a.control}: {'PASS' if result['passed'] else 'FAIL'} {result['verdict']}")
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(result, fh, indent=2)
    sys.exit(0 if (not a.control or result["passed"]) else 1)


if __name__ == "__main__":
    main()
