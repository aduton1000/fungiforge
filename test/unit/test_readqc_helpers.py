"""Unit tests for the W2.2 read-QC helpers: bin/kmer_profile.py (GenomeScope2 parsing, ploidy
hint, coverage), bin/read_stats.py (fastp / NanoPlot parsing) and bin/read_triage.py."""
import json, os, subprocess, sys
import kmer_profile as kp  # noqa: E402  (sys.path set in conftest)
import read_stats as rs  # noqa: E402
import read_triage as rt  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

GS_P1 = """GenomeScope version 2.0
input file = kmc.hist
output directory = gs_p1
p = 1
k = 21

property                      min               max
Homozygous (a)                100%              100%
Genome Haploid Length         28,833,921 bp     28,894,742 bp
Genome Repeat Length          1,100,000 bp      1,120,000 bp
Genome Unique Length          27,700,000 bp     27,780,000 bp
Model Fit                     96.1412%          98.5103%
Read Error Rate               0.190513%         0.190513%
"""
GS_P2 = """GenomeScope version 2.0
p = 2
k = 21

property                      min               max
Homozygous (aa)               99.4536%          99.5257%
Heterozygous (ab)             0.474312%         0.546358%
Genome Haploid Length         28,700,000 bp     28,800,000 bp
Genome Repeat Length          1,000,000 bp      1,050,000 bp
Genome Unique Length          27,700,000 bp     27,750,000 bp
Model Fit                     94.0%             97.9%
Read Error Rate               0.2%              0.2%
"""
MODEL = """Formula: y ~ ...
Parameters:
         Estimate Std. Error t value Pr(>|t|)
d       1.234e-02  1.0e-03    12.3   <2e-16 ***
r       4.700e-03  1.0e-04    47.0   <2e-16 ***
kmercov 3.512e+01  1.0e-01   351.2   <2e-16 ***
bias    5.0e-01    1.0e-02    50.0   <2e-16 ***
length  2.886e+07  1.0e+04  2886.0   <2e-16 ***
"""


def gs_dir(tmp_path, name, summary, model=MODEL):
    d = tmp_path / name; d.mkdir()
    (d / "summary.txt").write_text(summary)
    if model:
        (d / "model.txt").write_text(model)
    return str(d)


# ---------------------------------------------------------------- kmer_profile
def test_parse_summary_and_model(tmp_path):
    d = gs_dir(tmp_path, "p1", GS_P1)
    s = kp.parse_summary(os.path.join(d, "summary.txt"))
    assert s["Genome Haploid Length"] == (28833921.0, 28894742.0)
    assert s["Model Fit"] == (96.1412, 98.5103) and s["Homozygous (a)"] == (100.0, 100.0)
    assert kp.parse_model_kcov(os.path.join(d, "model.txt")) == 35.12
    fit = kp.read_fit(d)
    assert fit["genome_size"] == 28864332 and fit["heterozygosity_pct"] is None and fit["converged"] and fit["kcov"] == 35.12
    assert kp.read_fit(str(tmp_path / "absent")) is None


def test_ploidy_rule():
    p1 = {"converged": True, "model_fit_pct": 96.1, "heterozygosity_pct": None}
    p2_low = {"converged": True, "model_fit_pct": 94.0, "heterozygosity_pct": 0.51}
    p2_het = {"converged": True, "model_fit_pct": 97.0, "heterozygosity_pct": 1.2}
    assert kp.choose_ploidy(p1, p2_low)[0] == "haploid"                   # p2 fits worse
    assert kp.choose_ploidy(p1, p2_het)[0] == "diploid"                   # p2 fits better and het >= 0.5 %
    assert kp.choose_ploidy(None, {"converged": True, "model_fit_pct": 90.0, "heterozygosity_pct": 0.1})[0] == "undetermined"
    assert kp.choose_ploidy(None, None)[0] == "undetermined"
    assert kp.choose_ploidy({"converged": False, "model_fit_pct": None, "heterozygosity_pct": None}, None)[0] == "undetermined"


def test_hist_totals(tmp_path):
    h = tmp_path / "h.txt"; h.write_text("1\t1000\n2\t500\n35\t20000\n70\t300\n")
    assert kp.hist_totals(str(h)) == (21800, 1000 + 1000 + 700000 + 21000)
    assert kp.hist_totals(str(h), min_count=2) == (20800, 1000 + 700000 + 21000)


def test_kmer_profile_cli_reports_haploid_genome_and_coverage(tmp_path):
    p1 = gs_dir(tmp_path, "p1", GS_P1); p2 = gs_dir(tmp_path, "p2", GS_P2)
    h = tmp_path / "h.txt"; h.write_text("1\t1000\n35\t20000\n")
    out = tmp_path / "k.json"
    r = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "kmer_profile.py"), "--sample", "S1", "--hist", str(h),
                        "--gs-p1", p1, "--gs-p2", p2, "--read-bases", "1443216600", "--mean-read-len", "150",
                        "--source", "illumina_trimmed", "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(out))
    assert d["stage"] == "kmer" and d["ploidy_hint"] == "haploid" and d["genome_size_est"] == 28864332
    assert d["kmer_coverage"] == 35.12 and d["coverage_from_kmers"] == round(35.12 * 150 / 130, 1)
    assert d["coverage_from_read_bases"] == round(1443216600 / 28864332, 1)
    assert d["fits"]["p2"]["heterozygosity_pct"] == round((0.474312 + 0.546358) / 2, 4)
    assert "28.86 Mb" in r.stdout


def test_kmer_profile_cli_without_any_fit_is_honest(tmp_path):
    out = tmp_path / "k.json"
    r = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "kmer_profile.py"), "--sample", "S1",
                        "--gs-p1", str(tmp_path / "none1"), "--gs-p2", str(tmp_path / "none2"), "--out", str(out)], capture_output=True, text=True)
    d = json.load(open(out))
    assert r.returncode == 0 and d["genome_size_est"] is None and d["ploidy_hint"] == "undetermined" and "no GenomeScope model" in d["note"]


# ---------------------------------------------------------------- read_stats
FASTP = {"summary": {"before_filtering": {"total_reads": 26000000, "total_bases": 3926000000},
                     "after_filtering": {"total_reads": 25400000, "total_bases": 3800000000, "q30_rate": 0.93, "gc_content": 0.48}},
         "duplication": {"rate": 0.12}, "adapter_cutting": {"adapter_trimmed_reads": 1200000}}
NANO = """General summary:
Mean read length:                  9,530.2
Mean read quality:                    14.9
Median read length:                6,101.0
Median read quality:                  15.3
Number of reads:                 401,000.0
Read length N50:                  16,872.0
STDEV read length:                 8,700.0
Total bases:             3,821,000,000.0
Number, percentage and megabases of reads above quality cutoffs
>Q10:	392,000 (97.8%) 3760.1Mb
>Q15:	220,000 (54.9%) 2100.0Mb
"""


def test_fastp_and_nanoplot_parsing(tmp_path):
    fp = tmp_path / "fastp.json"; fp.write_text(json.dumps(FASTP))
    st = rs.fastp_stats(str(fp))
    assert st["reads_trimmed"] == 25400000 and st["bases_raw"] == 3926000000 and st["q30_rate_trimmed"] == 0.93
    assert round(st["mean_len_trimmed"], 1) == round(3800000000 / 25400000, 1) and st["duplication_rate"] == 0.12
    ns = tmp_path / "NanoStats.txt"; ns.write_text(NANO)
    n = rs.nanostats(str(ns))
    assert n["reads"] == 401000 and n["n50"] == 16872 and n["bases"] == 3821000000 and n["mean_q"] == 14.9
    assert n["reads_over_q10_pct"] == 97.8 and n["reads_over_q15_pct"] == 54.9
    assert rs.fastp_stats(str(tmp_path / "absent.json")) is None and rs.nanostats(None) is None


def test_read_stats_cli_merges_into_existing_stage_json(tmp_path):
    fp = tmp_path / "fastp.json"; fp.write_text(json.dumps(FASTP))
    ns = tmp_path / "filtNanoStats.txt"; ns.write_text(NANO)
    existing = tmp_path / "S1.readqc.json"; existing.write_text(json.dumps({"sample": "S1", "stage": "readqc", "ont_filtered": "x.fq.gz"}))
    r = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "read_stats.py"), "--sample", "S1", "--platform", "hybrid",
                        "--fastp", str(fp), "--nanostats-filt", str(ns), "--merge-into", str(existing), "--out", str(existing)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(existing))
    assert d["ont_filtered"] == "x.fq.gz" and d["platform"] == "hybrid"
    assert d["read_stats"]["illumina"]["reads_trimmed"] == 25400000 and d["read_stats"]["ont_filtered"]["n50"] == 16872
    assert d["total_bases"] == 3800000000 + 3821000000 and d["ont_n50"] == 16872 and d["read_stats"]["ont_raw"] is None


# ---------------------------------------------------------------- read_triage
def test_read_triage_uses_shared_verdict_rules(tmp_path):
    rep = os.path.join(ROOT, "test", "fixtures", "kraken2", "k2.report")
    out = tmp_path / "t.json"
    r = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "read_triage.py"), "--sample", "S1", "--report", rep,
                        "--n-reads", "1000", "--total-reads", "13000000", "--platform", "illumina", "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(out))
    assert d["stage"] == "triage" and d["verdict"] == "fungal" and d["composition_pct"]["fungi"] == 36.0
    assert d["top_taxon"] == "Aspergillus fumigatus (30.0%)" and d["subsample_reads"] == 1000 and d["total_reads"] == 13000000


def test_read_triage_missing_report_is_not_run(tmp_path):
    out = tmp_path / "t.json"
    r = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "read_triage.py"), "--sample", "S1", "--report", str(tmp_path / "none"),
                        "--n-reads", "0", "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 0 and json.load(open(out))["verdict"] == "not_run"
    assert rt.summarize_report(os.path.join(ROOT, "test", "fixtures", "kraken2", "k2.report"))["verdict"] == "fungal"
