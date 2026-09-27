"""Unit tests for bin/benchmark_assembly.py — reference benchmark from minimap2 PAF (W1.1).
The metric functions are tested on hand-built PAF records; the end-to-end test builds a
synthetic reference and a fragmented, mutated, partly inverted assembly and is skipped when
minimap2 is not on PATH (it is installed in CI)."""
import json, os, random, shutil, subprocess, sys
import pytest
import benchmark_assembly as ba  # noqa: E402  (sys.path set in conftest)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def paf(q, qlen, qs, qe, strand, t, tlen, ts, te, cs="", tp="P"):
    alen = te - ts
    return "\t".join(map(str, [q, qlen, qs, qe, strand, t, tlen, ts, te, alen, alen, 60, f"tp:A:{tp}", f"cs:Z:{cs}"]))


def test_count_cs_and_nga50():
    assert ba.count_cs(":100*at:50+ggg:20-tt:5") == (1, 1, 1, 3, 2)
    assert ba.count_cs("") == (0, 0, 0, 0, 0)
    assert ba.nga50([500, 300, 200], 1000) == 500 and ba.nga50([400, 400, 200], 1000) == 400 and ba.nga50([], 10) == 0


def test_benchmark_metrics_on_hand_built_alignments():
    asm = {"c1": 6000, "c2": 3000, "c3": 500}
    ref = {"chr1": 10000}
    recs = ba.parse_paf([
        paf("c1", 6000, 0, 3000, "+", "chr1", 10000, 0, 3000, cs=":1000*ac:1999"),          # collinear pair -> no breakpoint
        paf("c1", 6000, 3000, 6000, "+", "chr1", 10000, 3000, 6000, cs=":3000"),
        paf("c2", 3000, 0, 1500, "+", "chr1", 10000, 6000, 7500, cs=":1500"),               # c2 second half inverted -> 1 breakpoint
        paf("c2", 3000, 1500, 3000, "-", "chr1", 10000, 7500, 9000, cs=":1490+aaaaaaaaaa"),
        paf("c1", 6000, 0, 1000, "+", "chr1", 10000, 9000, 10000, cs=":1000", tp="S"),      # secondary: ignored
    ])
    res = ba.benchmark(asm, ref, recs, min_block=1000)
    assert res["assembly_len"] == 9500 and res["n_contigs"] == 3 and res["reference_len"] == 10000
    assert res["assembly_aligned_frac"] == round(9000 / 9500, 4) and res["reference_covered_frac"] == 0.9
    assert res["nga50"] == 3000 and res["breakpoints"] == 1 and res["contigs_with_breakpoints"] == 1
    assert res["mismatches"] == 1 and res["indels"] == 1 and res["inserted_bp"] == 10
    assert res["unaligned_contigs"] == 1 and res["unaligned_bp"] == 500
    assert res["per_contig"]["c2"]["breakpoints"] == 1 and res["per_contig"]["c1"]["blocks"] == 2


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not on PATH")
def test_end_to_end_synthetic_genome(tmp_path, helpers):
    rng = random.Random(7)
    ref = helpers["random_dna"](300000, seed=7)
    # assembly: three contigs; contig B carries 30 substitutions; contig C is inverted in its 2nd half; 3 kb missing between A and B
    A = ref[:100000]
    B = list(ref[103000:200000]); pos = rng.sample(range(1000, 96000), 30)
    for p in pos:
        B[p] = "A" if B[p] != "A" else "C"
    B = "".join(B)
    C = ref[200000:250000] + revcomp(ref[250000:300000])
    helpers["write_fasta"](tmp_path / "ref.fa", {"chr1": ref}, width=80)
    helpers["write_fasta"](tmp_path / "asm.fa", {"A": A, "B": B, "C": C, "junk": helpers["random_dna"](800, seed=99)}, width=80)
    out = tmp_path / "bench.json"
    r = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "benchmark_assembly.py"), "--assembly", str(tmp_path / "asm.fa"),
                        "--reference", str(tmp_path / "ref.fa"), "--out", str(out), "--sample", "synthetic"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    res = json.load(open(out))
    assert res["n_contigs"] == 4 and res["reference_len"] == 300000
    assert res["reference_covered_frac"] >= 0.98                         # 3 kb gap only
    assert res["breakpoints"] == 1 and res["per_contig"]["C"]["breakpoints"] == 1
    assert 25 <= res["mismatches"] <= 35                                 # the 30 planted substitutions
    assert res["unaligned_contigs"] == 1 and res["unaligned_bp"] == 800
    assert res["nga50"] == 100000 or res["nga50"] >= 97000
    assert "NGA50" in r.stdout
