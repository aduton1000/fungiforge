#!/usr/bin/env python3
"""K-mer profile of an isolate from KMC + GenomeScope2 outputs (W2.2, stage 01b).

  kmer_profile.py --sample ID --hist kmc.hist --gs-p1 DIR --gs-p2 DIR [--k 21]
                  [--read-bases N --mean-read-len L] --out ID.kmer.json

GenomeScope2 is run twice by the module (ploidy 1 and 2); this script parses both
`summary.txt` / `model.txt`, picks a ploidy hint, and reports genome size, heterozygosity,
k-mer coverage and read coverage. Anything that could not be fitted is reported as null with
a note — never invented. Pure stdlib.

Ploidy hint rule (conservative): diploid when the p=2 model converged, its fit is at least as
good as the p=1 fit, and its heterozygosity is >= --het-diploid (default 0.5 %); haploid when the
p=1 model converged; otherwise "undetermined".
"""
from __future__ import annotations
import argparse, json, os, re

NUM = re.compile(r"[-+]?\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?")


def _num(s):
    m = NUM.search(s or "")
    return float(m.group(0).replace(",", "")) if m else None


def parse_summary(path):
    """summary.txt -> {property: (min, max)} with numbers parsed (percent -> float, bp -> int-ish)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    header_seen = False
    for line in open(path):
        line = line.rstrip("\n")
        if line.lower().startswith("property"):
            header_seen = True
            continue
        if not header_seen or not line.strip():
            continue
        # property name = text before the first run of 2+ spaces; then min, max
        m = re.match(r"^(.*?)\s{2,}(\S+(?: bp)?)\s+(\S+(?: bp)?)\s*$", line)
        if not m:
            continue
        prop, lo, hi = m.group(1).strip(), m.group(2), m.group(3)
        out[prop] = (_num(lo), _num(hi))
    return out


def parse_model_kcov(path):
    """kmercov from model.txt (the nls parameter block); None when absent."""
    if not path or not os.path.exists(path):
        return None
    for line in open(path):
        if re.match(r"^\s*kmercov\b", line):
            v = _num(line.split()[1]) if len(line.split()) > 1 else None
            return v
    return None


def read_fit(gs_dir):
    """One GenomeScope run -> dict (or None when the run produced no summary)."""
    summ = parse_summary(os.path.join(gs_dir or "", "summary.txt"))
    if not summ:
        return None
    def g(name):
        for k, v in summ.items():
            if k.lower().startswith(name.lower()):
                return v
        return (None, None)
    hap = g("Genome Haploid Length"); het = g("Heterozygous"); fit = g("Model Fit"); err = g("Read Error Rate")
    rep = g("Genome Repeat Length"); uniq = g("Genome Unique Length")
    kcov = parse_model_kcov(os.path.join(gs_dir, "model.txt"))
    return {"genome_size_min": hap[0], "genome_size_max": hap[1],
            "genome_size": round((hap[0] + hap[1]) / 2) if hap[0] is not None and hap[1] is not None else None,
            "heterozygosity_pct": round((het[0] + het[1]) / 2, 4) if het[0] is not None and het[1] is not None else None,
            "model_fit_pct": fit[0], "read_error_pct": err[0],
            "repeat_len": rep[1], "unique_len": uniq[1], "kcov": kcov,
            "converged": fit[0] is not None and hap[0] is not None}


def choose_ploidy(p1, p2, het_diploid=0.5):
    """Return (ploidy_hint, chosen_fit_key, note)."""
    ok1 = bool(p1 and p1["converged"]); ok2 = bool(p2 and p2["converged"])
    if ok2 and (not ok1 or (p2["model_fit_pct"] or 0) >= (p1["model_fit_pct"] or 0)) \
            and (p2["heterozygosity_pct"] or 0) >= het_diploid:
        return "diploid", "p2", f"p=2 model fits at least as well and heterozygosity {p2['heterozygosity_pct']} % >= {het_diploid} %"
    if ok1:
        return "haploid", "p1", "p=1 model converged" + ("" if not ok2 else f"; p=2 heterozygosity {p2['heterozygosity_pct']} % below the diploid threshold")
    if ok2:
        return "undetermined", "p2", "only the p=2 model converged but heterozygosity is below the diploid threshold"
    return "undetermined", None, "no GenomeScope model converged (low coverage, very noisy reads, or contamination)"


def hist_totals(hist_path, min_count=1):
    """Total distinct k-mers and total k-mer instances from a KMC histogram (count\tfreq)."""
    distinct = total = 0
    if hist_path and os.path.exists(hist_path):
        for line in open(hist_path):
            f = line.split()
            if len(f) >= 2 and f[0].isdigit():
                c, n = int(f[0]), int(f[1])
                if c >= min_count:
                    distinct += n; total += c * n
    return distinct, total


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--hist", default=None)
    ap.add_argument("--gs-p1", default=None); ap.add_argument("--gs-p2", default=None)
    ap.add_argument("--k", type=int, default=21); ap.add_argument("--het-diploid", type=float, default=0.5)
    ap.add_argument("--read-bases", type=float, default=None, help="total bases in the reads profiled")
    ap.add_argument("--mean-read-len", type=float, default=None)
    ap.add_argument("--source", default="reads", help="which reads were profiled (illumina_trimmed | ont_filtered)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    p1, p2 = read_fit(a.gs_p1), read_fit(a.gs_p2)
    ploidy, key, note = choose_ploidy(p1, p2, a.het_diploid)
    chosen = {"p1": p1, "p2": p2}.get(key) if key else None
    distinct, total = hist_totals(a.hist)
    genome_size = chosen["genome_size"] if chosen else None
    kcov = chosen["kcov"] if chosen else None
    cov_kmer = None
    if kcov and a.mean_read_len and a.mean_read_len > a.k:
        cov_kmer = round(kcov * a.mean_read_len / (a.mean_read_len - a.k + 1), 1)
    cov_raw = round(a.read_bases / genome_size, 1) if (a.read_bases and genome_size) else None
    out = {"sample": a.sample, "stage": "kmer", "k": a.k, "source": a.source,
           "genome_size_est": genome_size,
           "heterozygosity_pct": chosen["heterozygosity_pct"] if chosen else None,
           "ploidy_hint": ploidy, "ploidy_basis": note,
           "kmer_coverage": kcov, "coverage_from_kmers": cov_kmer, "coverage_from_read_bases": cov_raw,
           "read_bases": a.read_bases, "mean_read_len": a.mean_read_len,
           "hist_distinct_kmers": distinct, "hist_total_kmers": total,
           "fits": {"p1": p1, "p2": p2}}
    if not chosen:
        out["note"] = note
    json.dump(out, open(a.out, "w"), indent=2)
    gs = f"{genome_size/1e6:.2f} Mb" if genome_size else "NA"
    print(f"[kmer_profile] {a.sample}: genome ~{gs}, het {out['heterozygosity_pct']} %, ploidy hint {ploidy}, "
          f"k-mer cov {kcov}, coverage ~{cov_kmer or cov_raw} ({a.source})")


if __name__ == "__main__":
    main()
