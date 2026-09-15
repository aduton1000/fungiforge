"""Unit tests for bin/kraken_taxonomy.py — report parsing, lineage grouping, verdict rules,
and the classify-contigs / summarize sub-commands (fixture: test/fixtures/kraken2/k2.report)."""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
BIN = os.path.join(ROOT, "bin")
REPORT = os.path.join(ROOT, "test", "fixtures", "kraken2", "k2.report")
sys.path.insert(0, BIN)
import kraken_taxonomy as kt  # noqa: E402


def test_parse_report_reconstructs_parents_from_indentation():
    nodes = kt.parse_report(REPORT)
    assert nodes[746128]["parent"] == 5052 and nodes[5052]["parent"] == 4751
    assert nodes[4751]["parent"] == 2759 and nodes[2759]["parent"] == 131567
    assert nodes[562]["parent"] == 2 and nodes[2]["parent"] == 131567
    assert nodes[0]["parent"] is None and nodes[1]["parent"] is None
    assert nodes[746128]["rank"] == "S" and nodes[746128]["clade_pct"] == 30.0 and nodes[746128]["taxon_n"] == 300


def test_group_of_assigns_domains_human_and_fungi():
    nodes = kt.parse_report(REPORT)
    assert kt.group_of(nodes, 746128) == "fungi"
    assert kt.group_of(nodes, 5052) == "fungi"
    assert kt.group_of(nodes, 9606) == "human"
    assert kt.group_of(nodes, 562) == "bacteria"
    assert kt.group_of(nodes, 33208) == "other_eukaryota"
    assert kt.group_of(nodes, 131567) == "unresolved"
    assert kt.group_of(nodes, 0) == "unclassified" and kt.group_of(nodes, 999999) == "unclassified"


def test_verdict_rules():
    assert kt.verdict({"human": 0.6}) == "human"
    assert kt.verdict({"bacteria": 0.4, "viruses": 0.1}) == "non_fungal"
    assert kt.verdict({"fungi": 0.15, "bacteria": 0.1, "unclassified": 0.75}) == "fungal"
    assert kt.verdict({"unclassified": 0.9, "fungi": 0.05, "bacteria": 0.05}) == "likely_fungal"
    assert kt.verdict({"unclassified": 0.6, "bacteria": 0.15, "fungi": 0.05}) == "mixed"
    assert kt.verdict({"fungi": 0.3, "bacteria": 0.3}) == "mixed"


def test_top_species_sorted_by_clade_pct_with_group():
    nodes = kt.parse_report(REPORT)
    top = kt.top_species(nodes)
    assert [t["name"] for t in top] == ["Aspergillus fumigatus", "Escherichia coli", "Aspergillus flavus"]
    assert top[0]["group"] == "fungi" and top[1]["group"] == "bacteria"


def run(*args):
    return subprocess.run([sys.executable, os.path.join(BIN, "kraken_taxonomy.py"), *args], capture_output=True, text=True)


def test_summarize_json_reports_composition_and_verdict():
    r = run("summarize", REPORT, "--json")
    assert r.returncode == 0, r.stderr
    d = json.loads(r.stdout)
    assert d["total"] == 1000 and d["pct"]["fungi"] == 36.0 and d["pct"]["bacteria"] == 10.0 and d["pct"]["human"] == 4.0
    assert d["verdict"] == "fungal"


def test_summarize_tsv_has_eleven_columns():
    r = run("summarize", REPORT)
    f = r.stdout.strip().split("\t")
    assert len(f) == 11 and f[0] == "1000" and f[-1] == "fungal" and f[8] == "Aspergillus fumigatus"


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n{s}\n")


def test_classify_contigs_drops_bacterial_and_human_contigs(tmp_path):
    asm = tmp_path / "asm.fa"
    write_fasta(asm, {"c_fungal": "A" * 900, "c_bact": "C" * 50, "c_human": "G" * 30, "c_unc": "T" * 20})
    k2out = tmp_path / "k2.out"
    k2out.write_text("C\tc_fungal\t746128\t900\tx\nC\tc_bact\t562\t50\tx\nC\tc_human\t9606\t30\tx\nU\tc_unc\t0\t20\tx\n")
    r = run("classify-contigs", REPORT, str(k2out), str(asm), "--sample", "S1",
            "--out-fasta", str(tmp_path / "kept.fa"), "--out-json", str(tmp_path / "d.json"), "--mito", "S1.mito.fasta")
    assert r.returncode == 0, r.stderr
    d = json.load(open(tmp_path / "d.json"))
    kept = [l[1:].strip() for l in open(tmp_path / "kept.fa") if l.startswith(">")]
    assert kept == ["c_fungal", "c_unc"]
    assert d["verdict"] == "fungal" and d["filtered"] is True and d["n_contigs_kept"] == 2
    assert d["bp_in"] == 1000 and d["bp_kept"] == 920 and d["dropped_bp_pct"] == 8.0
    assert d["composition_bp_pct"] == {"fungi": 90.0, "bacteria": 5.0, "human": 3.0, "unclassified": 2.0}
    assert d["stage"] == "decontam" and d["sample"] == "S1" and d["mito"] == "S1.mito.fasta"


def test_classify_contigs_keeps_a_bacterial_isolate_whole(tmp_path):
    asm = tmp_path / "asm.fa"
    write_fasta(asm, {"c1": "A" * 800, "c2": "C" * 200})
    k2out = tmp_path / "k2.out"; k2out.write_text("C\tc1\t562\t800\tx\nC\tc2\t746128\t200\tx\n")
    r = run("classify-contigs", REPORT, str(k2out), str(asm), "--out-fasta", str(tmp_path / "kept.fa"), "--out-json", str(tmp_path / "d.json"))
    assert r.returncode == 0, r.stderr
    d = json.load(open(tmp_path / "d.json"))
    assert d["verdict"] == "non_fungal" and d["filtered"] is False and d["n_contigs_kept"] == 2
    assert "kept whole" in d["note"]


def test_classify_contigs_keep_all_flag(tmp_path):
    asm = tmp_path / "asm.fa"; write_fasta(asm, {"c1": "A" * 900, "c2": "C" * 100})
    k2out = tmp_path / "k2.out"; k2out.write_text("C\tc1\t746128\t900\tx\nC\tc2\t562\t100\tx\n")
    r = run("classify-contigs", REPORT, str(k2out), str(asm), "--keep-all", "--out-fasta", str(tmp_path / "k.fa"), "--out-json", str(tmp_path / "d.json"))
    d = json.load(open(tmp_path / "d.json"))
    assert r.returncode == 0 and d["n_contigs_kept"] == 2 and d["filtered"] is False and d["note"].startswith("bacterial/")


def test_classify_contigs_wraps_sequences_at_80(tmp_path):
    asm = tmp_path / "asm.fa"; write_fasta(asm, {"c1": "A" * 170})
    k2out = tmp_path / "k2.out"; k2out.write_text("C\tc1\t746128\t170\tx\n")
    run("classify-contigs", REPORT, str(k2out), str(asm), "--out-fasta", str(tmp_path / "k.fa"), "--out-json", str(tmp_path / "d.json"))
    lines = open(tmp_path / "k.fa").read().splitlines()
    assert lines == [">c1", "A" * 80, "A" * 80, "A" * 10]
