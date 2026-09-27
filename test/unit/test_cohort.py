"""Unit tests for bin/busco_singlecopy.py and bin/cohort_phylo.py (W3.1): BUSCO/compleasm layouts,
ANI clustering, paftools call parsing, core intersection, SNP distances, clonal groups, trimming,
concatenation, gene selection, and the CLI without external tools."""
import json, os, subprocess, sys
import busco_singlecopy as bs  # noqa: E402  (sys.path set in conftest)
import cohort_phylo as cp  # noqa: E402


def test_busco_singlecopy_both_layouts(tmp_path, helpers):
    b = tmp_path / "busco" / "run_fungi_odb10" / "busco_sequences" / "single_copy_busco_sequences"; b.mkdir(parents=True)
    (b / "1234at4751.faa").write_text(">c1:100-400\nMAAA*\n"); (b / "5678at4751.faa").write_text(">c2:1-90\nMCCC\n")
    assert bs.from_busco(str(tmp_path / "busco")) == {"1234at4751": "MAAA*", "5678at4751": "MCCC"}
    c = tmp_path / "compleasm" / "fungi_odb10"; c.mkdir(parents=True)
    (c / "full_table.tsv").write_text("# Gene\tStatus\tSequence\n1234at4751\tSingle\tc1\n9999at4751\tDuplicated\tc3\n5678at4751\tMissing\t\n")
    (c / "translated_protein.fasta").write_text(">1234at4751_1\nMDDD\n>9999at4751_1\nMEEE\n>9999at4751_2\nMFFF\n")
    assert bs.from_compleasm(str(tmp_path / "compleasm")) == {"1234at4751": "MDDD"}
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "busco_singlecopy.py"), "--sample", "S1", "--busco-dir", str(tmp_path / "busco"),
                        "--compleasm-dir", str(tmp_path / "compleasm"), "--out", str(tmp_path / "sc.faa")], capture_output=True, text=True)
    assert r.returncode == 0 and (tmp_path / "sc.faa").read_text() == ">1234at4751|S1\nMDDD\n"          # compleasm wins when present
    r2 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "busco_singlecopy.py"), "--sample", "S1", "--busco-dir", str(tmp_path / "busco"),
                         "--compleasm-dir", str(tmp_path / "none"), "--out", str(tmp_path / "sc2.faa")], capture_output=True, text=True)
    assert r2.returncode == 0 and (tmp_path / "sc2.faa").read_text() == ">1234at4751|S1\nMAAA\n>5678at4751|S1\nMCCC\n"
    r3 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "busco_singlecopy.py"), "--sample", "S1", "--busco-dir", str(tmp_path / "x"), "--compleasm-dir", str(tmp_path / "y"), "--out", str(tmp_path / "sc3.faa")], capture_output=True, text=True)
    assert r3.returncode == 0 and (tmp_path / "sc3.faa").read_text() == "" and "no single-copy" in r3.stderr


def test_sample_of_and_ani_clusters():
    assert cp.sample_of("/x/asm/S1.nuclear.fasta") == "S1" and cp.sample_of("busco/S2.busco_sc.faa") == "S2" and cp.sample_of("GCF_1.fna") == "GCF_1"
    txt = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n" \
          "A.nuclear.fasta\tB.nuclear.fasta\t99.5\t90\t91\tA\tB\nA.nuclear.fasta\tC.nuclear.fasta\t88.0\t70\t70\tA\tC\nB.nuclear.fasta\tC.nuclear.fasta\t87.9\t70\t70\tB\tC\n" \
          "A.nuclear.fasta\tD.nuclear.fasta\t99.0\t3\t2\tA\tD\n"
    e = cp.parse_skani_edges(txt)
    assert e[("A", "B")] == (99.5, 91.0) and e[("A", "D")] == (99.0, 3.0)
    assert cp.ani_clusters(["A", "B", "C", "D"], e) == [["A", "B"], ["C"], ["D"]]        # D: ANI high but aligned fraction too low
    assert cp.ani_clusters(["A", "B", "C"], e, threshold=85) == [["A", "B", "C"]]
    assert cp.ani_clusters(["X"], {}) == [["X"]]


def test_paftools_call_parsing_core_and_snp_distances():
    call_b = "R\tchr1\t0\t60000\nV\tchr1\t10000\t10001\t1\t60\tt\ta\tB\t10000\t10001\t+\nV\tchr1\t45000\t45010\t1\t60\tagcgcaactc\t-\tB\t45000\t45000\t+\n"
    call_c = "R\tchr1\t20000\t60000\nV\tchr1\t30000\t30001\t1\t60\ta\tt\tC\t30000\t30001\t+\nV\tchr1\t10000\t10001\t1\t60\tt\tg\tC\t9000\t9001\t+\n"
    rb, sb, ib = cp.parse_paftools_call(call_b)
    assert rb == {"chr1": [(0, 60000)]} and sb == {("chr1", 10000): "A"} and ib == 1
    rc, sc, _ = cp.parse_paftools_call(call_c)
    core = cp.intersect_regions([rb, rc])
    assert core == {"chr1": [(20000, 60000)]}
    core_bp, d = cp.snp_distances(["A", "B", "C"], {"B": (rb, sb), "C": (rc, sc)}, "A")
    # B's SNV at 10000 and C's at 10000 lie outside the core; C's SNV at 30000 is inside
    assert core_bp == 40000 and d["A"]["B"] == 0 and d["A"]["C"] == 1 and d["B"]["C"] == 1 and d["B"]["B"] == 0
    groups = cp.single_linkage(["A", "B", "C"], lambda x, y: d[x][y] <= 0)
    assert groups == [["A", "B"], ["C"]]
    assert cp.choose_reference(["A", "B"], {"A": {"busco_complete": 98.0, "n50": 1e6}, "B": {"busco_complete": 99.5, "n50": 5e5}}) == "B"
    assert cp.choose_reference(["A", "B"], {}) == "A"


def test_trim_concatenate_and_gene_selection():
    aln = {"S1": "MA-KL", "S2": "MAAKL", "S3": "M--KL"}
    assert cp.trim_columns(aln) == {"S1": "MAKL", "S2": "MAKL", "S3": "M-KL"}        # column 3: 2/3 gaps dropped; column 2: 1/3 kept
    sm, parts = cp.concatenate([{"S1": "AAA", "S2": "AAC"}, {"S1": "GG"}], ["S1", "S2"])
    assert sm == {"S1": "AAAGG", "S2": "AAC--"} and parts == [("gene1", 1, 3), ("gene2", 4, 5)]
    sets = {"S1": {"g1": "M", "g2": "M", "g3": "M"}, "S2": {"g1": "M", "g2": "M"}, "S3": {"g1": "M"}, "S4": {"g1": "M", "g2": "M"}}
    assert cp.select_genes(sets, ["S1", "S2", "S3", "S4"], min_frac=0.75) == ["g1", "g2"]
    assert cp.select_genes(sets, ["S1", "S2", "S3", "S4"], min_frac=1.0) == ["g1"]


def test_cli_without_tools_writes_clusters_and_json(tmp_path, helpers):
    asm = [helpers["write_fasta"](tmp_path / f"{s}.nuclear.fasta", {"c1": helpers["random_dna"](2000, seed=i)}) for i, s in enumerate(["S1", "S2"])]
    sp = []
    for s in ("S1", "S2"):
        p = tmp_path / f"{s}.species.txt"; p.write_text("Aspergillus fumigatus\thigh\n"); sp.append(str(p))
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "cohort_phylo.py"), "--outdir", str(tmp_path / "cohort"), "--assemblies", *asm, "--species", *sp],
                       capture_output=True, text=True, env={**os.environ, "PATH": "/nonexistent"})
    assert r.returncode == 0, r.stderr
    d = json.load(open(tmp_path / "cohort" / "cohort.json"))
    assert d["n_isolates"] == 2 and d["species_clusters"] == [{"cluster": "C1", "members": ["S1"], "species": ["Aspergillus fumigatus"]}, {"cluster": "C2", "members": ["S2"], "species": ["Aspergillus fumigatus"]}]
    assert d["phylogenomics"]["tree"] is None and "needs >= 4" in d["phylogenomics"]["note"] and d["clusters_snp"] == []
    assert (tmp_path / "cohort" / "species_clusters.tsv").read_text().splitlines()[0] == "cluster\tsample\tspecies"
    assert (tmp_path / "cohort" / "clonal_groups.tsv").exists()
