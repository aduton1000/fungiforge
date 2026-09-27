#!/usr/bin/env python3
"""Read-level genotyping of resistance hotspots (W2.4, stage 09).

The assembly-level caller (af_resistance.py) reads residues from the annotated proteins; this
module confirms them from the reads themselves: every read that covers a hotspot codon of the
isolate's own gene (CDS coordinates from the Funannotate GenBank) is translated at that codon,
giving per-residue allele counts, allele frequency and zygosity (a diploid Candida can carry a
resistance allele heterozygously, which an assembly collapses), and the cyp51A promoter tandem
repeat is checked as an insertion in the reads spanning the TR site. Copy-number signal for
target/efflux genes comes from locus depth over genome depth.

Library use (from af_resistance.py) or CLI:
  read_genotype.py --bam reads.bam --gbk isolate.gbk --protein FUN_000123-T1 --protein-fasta prot.faa
                   --residues 98,121,289 [--out geno.json]

`--sam <file>` replaces the BAM (tests): the whole SAM is parsed once instead of samtools view.
Reads must be mapped to the sequences of the GenBank records (bin/gbk_to_fasta.py) so that
GenBank coordinates are the BAM coordinates.
"""
from __future__ import annotations
import argparse, collections, json, os, re, subprocess, sys

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODON = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_AA):
    CODON[_B[_i // 16] + _B[(_i // 4) % 4] + _B[_i % 4]] = _a


def revcomp(s):
    return s.translate(COMP)[::-1]


def complement(s):
    return s.translate(COMP)


def translate_codon(c):
    c = c.upper()
    if "-" in c:
        return "del"
    return CODON.get(c, "X")


def parse_cigar(cigar):
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def walk(pos, ops, seq):
    """(ref_pos_1based, op, length, payload) per CIGAR element (I: position it follows)."""
    r, q = pos, 0
    for n, op in ops:
        if op in ("M", "=", "X"):
            for i in range(n):
                yield r + i, "M", 1, seq[q + i] if q + i < len(seq) else "N"
            r += n; q += n
        elif op == "I":
            yield r - 1, "I", n, seq[q:q + n]; q += n
        elif op == "D":
            yield r, "D", n, ""; r += n
        elif op == "N":
            r += n
        elif op == "S":
            q += n


# ---------- CDS coordinates from the GenBank ----------
def cds_codon_coords(gbk_path, protein_name=None, protein_seq=None):
    """Locate the CDS of `protein_name` (Funannotate: FUN_000123-T1 -> locus_tag FUN_000123 /
    protein_id ncbi:FUN_000123-T1) or, more robustly, the CDS whose /translation equals
    `protein_seq`. Returns {contig, strand, codons: [(p1, p2, p3), ...] (1-based genomic,
    transcript order), n_residues} or None."""
    try:
        from Bio import SeqIO
    except Exception:  # noqa: BLE001
        return None
    want = (protein_seq or "").rstrip("*").upper()
    base = re.sub(r"-T\d+$", "", protein_name or "")
    by_name = None
    for rec in SeqIO.parse(gbk_path, "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            names = " ".join(q.get("locus_tag", []) + q.get("protein_id", []) + q.get("transcript_id", []) + q.get("gene", []))
            tr = (q.get("translation", [""])[0] or "").rstrip("*").upper()
            hit = (rec.id, feat)
            if want and tr == want:
                return _codons(hit)
            if base and by_name is None and re.search(r"(^|[^A-Za-z0-9])" + re.escape(base) + r"([^A-Za-z0-9]|$)", names):
                by_name = hit
    return _codons(by_name) if by_name else None


def _codons(hit):
    rec_id, feat = hit
    strand = 1 if feat.location.strand in (1, None) else -1
    positions = []
    parts = sorted(feat.location.parts, key=lambda p: int(p.start))
    for p in parts:
        positions.extend(range(int(p.start) + 1, int(p.end) + 1))
    if strand == -1:
        positions.reverse()
    try:
        offset = int(feat.qualifiers.get("codon_start", ["1"])[0]) - 1
    except ValueError:
        offset = 0
    positions = positions[offset:]
    codons = [tuple(positions[i:i + 3]) for i in range(0, len(positions) - 2, 3)]
    return {"contig": rec_id, "strand": strand, "codons": codons, "n_residues": len(codons)}


# ---------- reads ----------
def sam_from_bam(bam, contig, start, end):
    cmd = ["samtools", "view", bam, f"{contig}:{start}-{end}"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"samtools view failed: {p.stderr.strip()[:200]}")
    return p.stdout.splitlines()


def _reads(sam_lines, contig, min_mapq=1):
    for line in sam_lines:
        if not line or line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11 or f[2] != contig:
            continue
        flag, pos, mapq, cigar, seq = int(f[1]), int(f[3]), int(f[4]), f[5], f[9]
        if flag & 0x904 or cigar == "*" or mapq < min_mapq:      # unmapped, secondary, supplementary
            continue
        yield pos, parse_cigar(cigar), seq


def codon_counts(sam_lines, contig, strand, codons_by_residue, min_mapq=1):
    """{residue: Counter(amino_acid)} from reads covering all three codon positions."""
    counts = {r: collections.Counter() for r in codons_by_residue}
    wanted = {p for trip in codons_by_residue.values() for p in trip}
    lo, hi = (min(wanted), max(wanted)) if wanted else (0, 0)
    for pos, ops, seq in _reads(sam_lines, contig, min_mapq):
        ref_len = sum(n for n, op in ops if op in "MDN=X")
        if pos > hi or pos + ref_len - 1 < lo:
            continue
        bases = {}
        for rpos, op, n, payload in walk(pos, ops, seq):
            if op == "M" and rpos in wanted:
                bases[rpos] = payload.upper()
            elif op == "D":
                for i in range(n):
                    if rpos + i in wanted:
                        bases[rpos + i] = "-"
        for r, trip in codons_by_residue.items():
            if all(p in bases for p in trip):
                codon = "".join(bases[p] for p in trip)   # trip is already in transcript order (descending on -)
                if strand == -1:
                    codon = complement(codon)
                counts[r][translate_codon(codon)] += 1
    return counts


def indel_support(sam_lines, contig, window, min_len, max_len, margin=20, min_mapq=1):
    """Reads spanning `window` (1-based inclusive, with margin) and how many carry an insertion
    or a deletion of min_len..max_len bp inside it: the read-level tandem-repeat check against
    the assembly used as mapping reference (an extra copy in the reads shows as an insertion, a
    copy the assembly has but the reads lack as a deletion)."""
    w0, w1 = window
    out = {"spanning": 0, "with_insertion": 0, "with_deletion": 0, "insertion_lengths": collections.Counter(),
           "deletion_lengths": collections.Counter()}
    for pos, ops, seq in _reads(sam_lines, contig, min_mapq):
        ref_len = sum(n for n, op in ops if op in "MDN=X")
        start, end = pos, pos + ref_len - 1
        if not (start <= w0 - margin and end >= w1 + margin):
            continue
        out["spanning"] += 1
        ins = [n for rpos, op, n, _ in walk(pos, ops, seq) if op == "I" and w0 - 5 <= rpos <= w1 + 5 and min_len <= n <= max_len]
        dels = [n for rpos, op, n, _ in walk(pos, ops, seq) if op == "D" and rpos <= w1 + 5 and rpos + n - 1 >= w0 - 5 and min_len <= n <= max_len]
        if ins:
            out["with_insertion"] += 1; out["insertion_lengths"][max(ins)] += 1
        if dels:
            out["with_deletion"] += 1; out["deletion_lengths"][max(dels)] += 1
    out["insertion_lengths"] = dict(sorted(out["insertion_lengths"].items()))
    out["deletion_lengths"] = dict(sorted(out["deletion_lengths"].items()))
    n = out["spanning"]
    out["ins_frac"] = round(out["with_insertion"] / n, 3) if n else None
    out["del_frac"] = round(out["with_deletion"] / n, 3) if n else None
    return out


def insertion_support(sam_lines, contig, window, min_len, max_len, margin=20, min_mapq=1):
    """Backwards-compatible view of indel_support (insertions only; `frac` = ins_frac)."""
    o = indel_support(sam_lines, contig, window, min_len, max_len, margin, min_mapq)
    return {"spanning": o["spanning"], "with_insertion": o["with_insertion"], "insertion_lengths": o["insertion_lengths"], "frac": o["ins_frac"]}


def zygosity(counter, min_reads=10, hom_frac=0.8, het_min=0.2):
    """Allele summary for one residue: homozygous | heterozygous | mixed | insufficient."""
    depth = sum(counter.values())
    if depth < min_reads:
        return {"depth": depth, "alleles": {}, "major": None, "major_frac": None, "call": "insufficient"}
    alleles = {aa: round(n / depth, 3) for aa, n in counter.most_common()}
    major, frac = counter.most_common(1)[0][0], counter.most_common(1)[0][1] / depth
    if frac >= hom_frac:
        call = "homozygous"
    else:
        second = counter.most_common(2)[1][1] / depth if len(counter) > 1 else 0.0
        call = "heterozygous" if second >= het_min else "mixed"
    return {"depth": depth, "alleles": alleles, "major": major, "major_frac": round(frac, 3), "call": call}


# ---------- depth / copy number ----------
def genome_mean_depth(bam):
    """Length-weighted mean depth over all contigs from `samtools coverage`."""
    p = subprocess.run(["samtools", "coverage", bam], capture_output=True, text=True)
    if p.returncode != 0:
        return None
    tot_bp = tot_depth = 0.0
    for line in p.stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        try:
            length, meandepth = int(f[2]) - int(f[1]) + 1, float(f[6])
        except (ValueError, IndexError):
            continue
        tot_bp += length; tot_depth += length * meandepth
    return round(tot_depth / tot_bp, 2) if tot_bp else None


def locus_mean_depth(bam, contig, start, end):
    p = subprocess.run(["samtools", "depth", "-a", "-r", f"{contig}:{start}-{end}", bam], capture_output=True, text=True)
    if p.returncode != 0:
        return None
    vals = [int(l.split("\t")[2]) for l in p.stdout.splitlines() if l.count("\t") >= 2]
    return round(sum(vals) / len(vals), 2) if vals else None


def copy_number(bam, contig, start, end, genome_depth=None):
    g = genome_depth if genome_depth is not None else genome_mean_depth(bam)
    loc = locus_mean_depth(bam, contig, start, end)
    if not g or loc is None:
        return {"locus_depth": loc, "genome_depth": g, "ratio": None}
    return {"locus_depth": loc, "genome_depth": g, "ratio": round(loc / g, 2)}


# ---------- one gene ----------
def genotype_residues(reads_source, cds, residues, min_reads=10, min_mapq=1):
    """reads_source: (bam_path) or ('sam', lines). Returns {residue: zygosity dict + codon coords}."""
    trip = {r: cds["codons"][r - 1] for r in residues if 1 <= r <= cds["n_residues"]}
    if not trip:
        return {}
    lo = min(p for t in trip.values() for p in t); hi = max(p for t in trip.values() for p in t)
    if isinstance(reads_source, tuple) and reads_source[0] == "sam":
        lines = reads_source[1]
    else:
        lines = sam_from_bam(reads_source, cds["contig"], lo, hi)
    counts = codon_counts(lines, cds["contig"], cds["strand"], trip, min_mapq)
    return {r: dict(zygosity(counts[r], min_reads), codon=list(trip[r])) for r in trip}


def load_sam(path):
    with open(path) as fh:
        return [l.rstrip("\n") for l in fh]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam"); ap.add_argument("--sam")
    ap.add_argument("--gbk", required=True); ap.add_argument("--protein", required=True)
    ap.add_argument("--protein-fasta", default=None, help="FASTA holding the protein (matched by translation)")
    ap.add_argument("--residues", required=True, help="comma list of 1-based residues (reference numbering = query numbering here)")
    ap.add_argument("--min-reads", type=int, default=10); ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if not (a.bam or a.sam):
        sys.exit("need --bam or --sam")
    pseq = None
    if a.protein_fasta and os.path.exists(a.protein_fasta):
        name, buf = None, []
        for line in open(a.protein_fasta):
            if line.startswith(">"):
                if name == a.protein:
                    break
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        pseq = "".join(buf) if name == a.protein else None
    cds = cds_codon_coords(a.gbk, a.protein, pseq)
    if not cds:
        sys.exit(f"CDS for {a.protein} not found in {a.gbk}")
    src = ("sam", load_sam(a.sam)) if a.sam else a.bam
    res = genotype_residues(src, cds, [int(x) for x in a.residues.split(",") if x.strip()], a.min_reads)
    doc = {"protein": a.protein, "contig": cds["contig"], "strand": cds["strand"], "residues": {str(k): v for k, v in res.items()}}
    if a.out:
        json.dump(doc, open(a.out, "w"), indent=2)
    for r, v in res.items():
        print(f"{a.protein}\t{r}\t{v['call']}\tdepth={v['depth']}\t{v['alleles']}")


if __name__ == "__main__":
    main()
