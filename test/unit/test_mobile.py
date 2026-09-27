"""Unit tests for the W2.8 stage-10 helpers: repeatmasker_tbl.py (real CEA10 table), mycovirus_screen.py,
mito_heg.py (six-frame ORFs, genetic code 4) and mobile_merge.py."""
import json, os, subprocess, sys
import repeatmasker_tbl as rt  # noqa: E402  (sys.path set in conftest)
import mycovirus_screen as ms  # noqa: E402
import mito_heg as mh  # noqa: E402
import mobile_merge as mm  # noqa: E402

TBL = os.path.join(os.path.dirname(__file__), "..", "fixtures", "repeatmasker", "AfumCEA10.tbl")


def test_parse_real_repeatmasker_table():
    d = rt.parse_tbl(TBL)
    assert d["n_sequences"] == 43 and d["total_length"] == 29983935 and d["gc_pct"] == 49.26
    assert d["te_percent"] == 4.67 and d["bases_masked"] == 1398996
    assert d["classes"]["LTR elements"] == {"n": 712, "bp": 877221, "pct": 2.93} and d["classes"]["Gypsy/DIRS1"]["pct"] == 2.29
    assert d["classes"]["DNA transposons"]["n"] == 177 and d["classes"]["Total interspersed repeats"]["pct"] == 4.03
    s = d["summary"]
    assert (s["te_percent"], s["ltr_pct"], s["gypsy_pct"], s["copia_pct"], s["dna_transposon_pct"], s["interspersed_pct"]) == (4.67, 2.93, 2.29, 0.32, 0.83, 4.03)
    assert rt.parse_tbl(TBL + ".missing") is None


def test_repeatmasker_tbl_cli(tmp_path, helpers):
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "repeatmasker_tbl.py"), "--tbl", TBL, "--out", str(tmp_path / "t.json")], capture_output=True, text=True)
    assert r.returncode == 0 and "te_percent=4.67" in r.stdout and json.load(open(tmp_path / "t.json"))["summary"]["ltr_pct"] == 2.93


def hit(contig, subj, pid, alen, qs, qe, title, bits=200.0, ev=1e-30):
    return "\t".join([contig, subj, str(pid), str(alen), str(qs), str(qe), "1", str(alen), str(ev), str(bits), title])


def test_mycovirus_screen_classifies_and_merges_loci(tmp_path, helpers):
    p = tmp_path / "h.tsv"
    p.write_text("\n".join([
        hit("c1", "acc|GENBANK|A1|GENBANK|N1|RNA-dependent RNA polymerase [Aspergillus fumigatus partitivirus 1]", 55.0, 300, 1000, 1900, "acc|GENBANK|A1|GENBANK|N1|RNA-dependent RNA polymerase [Aspergillus fumigatus partitivirus 1]", bits=400),
        hit("c1", "acc|GENBANK|A2|GENBANK|N2|capsid protein [Aspergillus fumigatus partitivirus 1]", 48.0, 200, 2200, 2800, "acc|GENBANK|A2|GENBANK|N2|capsid protein [Aspergillus fumigatus partitivirus 1]", bits=250),   # merged (gap < 1 kb)
        hit("c1", "acc|GENBANK|A3|GENBANK|N3|polyprotein [Drosophila melanogaster gypsy retrotransposon]", 60.0, 900, 50000, 52700, "acc|GENBANK|A3|GENBANK|N3|polyprotein [Drosophila melanogaster gypsy retrotransposon]"),
        hit("c2", "acc|GENBANK|A4|GENBANK|N4|major capsid protein [Escherichia phage T4]", 90.0, 500, 100, 1600, "acc|GENBANK|A4|GENBANK|N4|major capsid protein [Escherichia phage T4]"),
        hit("mt", "acc|GENBANK|A5|GENBANK|N5|RNA-dependent RNA polymerase [Rhizoctonia solani mitovirus IF]", 40.0, 150, 100, 550, "acc|GENBANK|A5|GENBANK|N5|RNA-dependent RNA polymerase [Rhizoctonia solani mitovirus IF]"),
        hit("c3", "acc|GENBANK|A6|GENBANK|N6|hypothetical [Hypoxylon fendleri mitovirus 1]", 30.0, 150, 100, 550, "acc|GENBANK|A6|GENBANK|N6|hypothetical [Hypoxylon fendleri mitovirus 1]"),   # identity too low
    ]) + "\n")
    rows = ms.parse_hits(str(p)); assert len(rows) == 6
    L = ms.loci(rows)
    assert len(L) == 4 and L[0]["contig"] == "c1" and L[0]["n_hits"] == 2 and L[0]["category"] == "mycovirus" and L[0]["best"] == "Aspergillus fumigatus partitivirus 1"
    assert [x["category"] for x in L] == ["mycovirus", "retroelement_like", "other_viral", "mycovirus"]
    assert ms.classify_title("x [Saccharomyces cerevisiae virus L-A]") == "mycovirus" and ms.classify_title("Ty1 polyprotein [yeast]") == "retroelement_like"
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "mycovirus_screen.py"), "--sample", "S1", "--hits", str(p), "--mito-contigs", "mt", "--out", str(tmp_path / "m.json")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(tmp_path / "m.json"))
    assert d["n_mycovirus_eve_candidates"] == 2 and d["retroelement_like_loci"] == 1 and d["other_viral_loci"] == 1
    assert {c["compartment"] for c in d["mycovirus_candidates"]} == {"nuclear", "mitochondrial"} and "DNA-only" in d["note"]
    assert json.load(open(tmp_path / "m.json"))["families"]["Aspergillus fumigatus partitivirus 1"] == 1


def test_six_frame_orfs_code4_and_heg_scan(tmp_path, helpers):
    # a 120-aa ORF on the minus strand whose codons include TGA (Trp in code 4, a stop in code 1)
    prot = ("MKLW" * 30)
    rev = {"M": "ATG", "K": "AAA", "L": "TTA", "W": "TGA"}
    cds = "".join(rev[a] for a in prot)
    seq = helpers["random_dna"](500, seed=21) + cds.translate(str.maketrans("ACGT", "TGCA"))[::-1] + "TTA" + helpers["random_dna"](400, seed=22)
    orfs = mh.six_frame_orfs({"mt": seq}, min_len=100)
    minus = [o for o in orfs if o["strand"] == "-" and prot in o["seq"]]
    assert len(minus) == 1 and minus[0]["contig"] == "mt" and minus[0]["aa"] >= 120 and minus[0]["id"].startswith("mt_-1_")
    assert mh.translate("TGA") == "W" and mh.translate("TAA") == "*"
    dt = tmp_path / "heg.domtbl"
    oid = minus[0]["id"]
    dt.write_text(" ".join([oid, "-", "120", "LAGLIDADG_1", "PF00961.20", "80", "1e-20", "60.0", "0.1", "1", "1", "1e-18", "1e-20", "50.0", "0.1", "1", "80", "1", "80", "1", "80", "0.9", "-"]) + "\n")
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "mito_heg.py"), "--sample", "S1", "--mito", helpers["write_fasta"](tmp_path / "m.fa", {"mt": seq}),
                        "--domtbl", str(dt), "--out", str(tmp_path / "h.json")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(tmp_path / "h.json"))
    assert d["status"] == "ok" and d["n_heg_orfs"] == 1 and d["heg_orfs"][0]["families"] == ["LAGLIDADG_1"] and d["by_family"] == {"LAGLIDADG_1": 1}
    r2 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "mito_heg.py"), "--sample", "S1", "--mito", str(tmp_path / "m.fa"), "--out", str(tmp_path / "h2.json")], capture_output=True, text=True)
    assert r2.returncode == 0 and json.load(open(tmp_path / "h2.json"))["status"] == "skipped"


def test_mobile_merge_with_genomad_summaries(tmp_path, helpers):
    g = tmp_path / "genomad_out" / "all_contigs_summary"; g.mkdir(parents=True)
    (g / "all_contigs_virus_summary.tsv").write_text("seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n"
                                                    "c1|provirus_1_5000\t5000\tProvirus\t1-5000\t5\t11\t0.9\tNA\t2\t1.0\tViruses;Duplodnaviria;Caudoviricetes\n"
                                                    "c9\t30000\tNo terminal repeats\tNA\t30\t11\t0.95\tNA\t5\t2.0\tViruses;Riboviria\n")
    (g / "all_contigs_plasmid_summary.tsv").write_text("seq_name\tlength\n" "c4\t8000\n")
    s = mm.genomad_summary(str(tmp_path / "genomad_out"))
    assert s == {"n_virus": 2, "n_plasmid": 1, "n_provirus": 1, "virus_taxa": {"Caudoviricetes": 1, "Riboviria": 1}, "ran": True}
    assert mm.genomad_summary(str(tmp_path / "none"))["ran"] is False
    tbl = tmp_path / "tbl.json"; tbl.write_text(json.dumps(rt.parse_tbl(TBL)))
    te = tmp_path / "te.json"; te.write_text(json.dumps({"n_families": 35, "by_class": {"LTR/Gypsy": 10}}))
    myco = tmp_path / "myco.json"; myco.write_text(json.dumps({"sample": "S1", "n_mycovirus_eve_candidates": 2}))
    heg = tmp_path / "heg.json"; heg.write_text(json.dumps({"sample": "S1", "n_heg_orfs": 4, "status": "ok"}))
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "mobile_merge.py"), "--sample", "S1", "--te", str(te), "--tbl", str(tbl), "--genomad", str(tmp_path / "genomad_out"),
                        "--mycovirus", str(myco), "--heg", str(heg), "--out", str(tmp_path / "mob.json")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(tmp_path / "mob.json"))
    assert d["te_percent"] == 4.67 and d["te_landscape"]["ltr_pct"] == 2.93 and d["n_te_families"] == 35
    assert d["n_genomad_virus"] == 2 and d["n_genomad_plasmid"] == 1 and d["n_mycovirus"] == 2 and d["n_mito_heg"] == 4
    r2 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "mobile_merge.py"), "--sample", "S1", "--out", str(tmp_path / "mob2.json")], capture_output=True, text=True)
    d2 = json.load(open(tmp_path / "mob2.json"))
    assert r2.returncode == 0 and d2["te_percent"] is None and d2["n_genomad_virus"] is None and d2["n_mycovirus"] is None
