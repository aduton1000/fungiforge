"""Unit tests for bin/id_classify.py — UNITE best-hit parsing, ITS thresholds, secondary-locus parsing (ties),
MLST parsing, the W2.3 concordance rule and the species-aware BUSCO lineage choice."""
import json, os, subprocess, sys
import id_classify as idc  # noqa: E402  (sys.path set in conftest)

UNITE = "Aspergillus_fumigatus|KC411932|SH1136727.09FU|refs|k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Aspergillus;s__Aspergillus_fumigatus"
UNITE_SP = "Aspergillus_sp|XX000001|SH0000001.09FU|reps|k__Fungi;p__Ascomycota;g__Aspergillus;s__Aspergillus_sp"


def b6(tmp_path, rows):
    p = tmp_path / "unite.b6"
    p.write_text("".join(f"its1\t{t}\t{pid}\t500\t2\t0\t1\t500\t1\t500\t0\t900\n" for t, pid in rows))
    return str(p)


def test_parse_unite_b6_takes_best_identity_and_splits_taxonomy(tmp_path):
    hit = idc.parse_unite_b6(b6(tmp_path, [(UNITE_SP, 97.0), (UNITE, 99.6), (UNITE_SP, 95.0)]))
    assert hit == {"pident": 99.6, "species": "Aspergillus fumigatus", "genus": "Aspergillus", "sh": "SH1136727.09FU"}


def test_parse_unite_b6_drops_sp_and_unidentified_species(tmp_path):
    hit = idc.parse_unite_b6(b6(tmp_path, [(UNITE_SP, 99.0)]))
    assert hit["species"] is None and hit["genus"] == "Aspergillus"
    assert idc.parse_unite_b6(None) is None and idc.parse_unite_b6(str(tmp_path / "none.b6")) is None


def test_top_gather_extracts_binomial(tmp_path):
    g = tmp_path / "g.csv"
    g.write_text("intersect_bp,f_orig_query,name\n1000,0.9,\"GCF_000002655.1 Aspergillus fumigatus Af293 strain=Af293\"\n")
    assert idc.top_gather(str(g))["species"] == "Aspergillus fumigatus"
    g.write_text("intersect_bp,name\n")
    assert idc.top_gather(str(g)) is None


LINEAGE_MAP = os.path.join(os.path.dirname(__file__), "..", "..", "fungiforge", "resources", "busco_lineages.tsv")


def locus_b6(tmp_path, locus, rows):
    """rows: (accession, pident, title[, alen, slen])"""
    p = tmp_path / f"{locus}.b6"
    lines = []
    for r in rows:
        acc, pid, title = r[:3]
        alen, slen = (r[3], r[4]) if len(r) > 4 else (500, 520)
        lines.append(f"q|{locus}\t{acc}\t{pid}\t{alen}\t1\t{alen}\t1\t{alen}\t1e-50\t{alen*1.8:.0f}\t900\t{slen}\t{title}\n")
    p.write_text("".join(lines))
    return str(p)


def run(tmp_path, h, its_text=">its1\nACGT\n", unite=None, gather=None, loci=None, mlst=None, lineage_map=None):
    its = tmp_path / "its.fa"; its.write_text(its_text)
    cmd = [sys.executable, os.path.join(h["BIN"], "id_classify.py"), "--sample", "S1", "--its", str(its),
           "--out-species", str(tmp_path / "sp.txt"), "--out-json", str(tmp_path / "id.json"),
           "--out-lineage", str(tmp_path / "lineage.txt")]
    if unite: cmd += ["--unite-b6", unite]
    if gather: cmd += ["--gather", gather]
    for name, path in (loci or {}).items():
        cmd += ["--locus-b6", f"{name}={path}"]
    if mlst: cmd += ["--mlst", mlst]
    if lineage_map: cmd += ["--lineage-map", lineage_map]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.load(open(tmp_path / "id.json")), open(tmp_path / "sp.txt").read().strip().split("\t")


def test_its_alone_is_medium_with_an_unavailable_flag(tmp_path, helpers):
    """W2.3: a species call needs a second line of evidence; ITS on its own is medium, not high."""
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 98.5)]))
    assert sp == ["Aspergillus fumigatus", "medium"] and doc["method"] == "ITS/UNITE(no secondary loci)" and doc["its_identity"] == 98.5
    assert doc["flags"] == ["secondary_loci_unavailable"] and doc["its_present"] is True and doc["stage"] == "identify"
    assert doc["busco_lineage"] == "fungi_odb10" and doc["busco_lineage_basis"] == "no map"


def test_its_plus_one_agreeing_locus_is_high(tmp_path, helpers):
    cam = locus_b6(tmp_path, "CaM", [("MH1.1", 99.8, "Aspergillus fumigatus strain CBS 133.61 calmodulin (CaM) gene"),
                                      ("MH2.1", 97.1, "Aspergillus lentulus strain CBS 117885 calmodulin (CaM) gene")])
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 99.6)]), loci={"CaM": cam}, lineage_map=LINEAGE_MAP)
    assert sp == ["Aspergillus fumigatus", "high"] and doc["method"] == "ITS+CaM(concordant)" and doc["flags"] == []
    assert doc["loci"]["CaM"]["level"] == "species" and doc["loci"]["CaM"]["pident"] == 99.8
    assert doc["concordance"]["secondary_agree"] == ["CaM"]
    assert doc["busco_lineage"] == "eurotiales_odb10" and doc["busco_lineage_basis"] == "genus:Aspergillus"
    assert open(tmp_path / "lineage.txt").read().strip() == "eurotiales_odb10"


def test_discordant_secondary_locus_drops_to_genus_with_flag(tmp_path, helpers):
    bena = locus_b6(tmp_path, "BenA", [("KX1.1", 99.9, "Aspergillus lentulus strain CBS 117885 beta-tubulin (benA) gene")])
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 99.0)]), loci={"BenA": bena})
    assert sp == ["Aspergillus sp.", "low"] and doc["method"].endswith("(discordant)")
    assert doc["flags"] == ["discordant:BenA=Aspergillus lentulus"]
    assert doc["concordance"]["secondary_disagree"] == ["BenA=Aspergillus lentulus"]


def test_tie_between_sibling_species_never_confirms_but_is_reported(tmp_path, helpers):
    """A. flavus / A. oryzae are identical on CaM: the locus is a tie, ITS stays medium with the candidates listed."""
    flavus = "Aspergillus_flavus|X|SH1.09FU|refs|k__Fungi;g__Aspergillus;s__Aspergillus_flavus"
    cam = locus_b6(tmp_path, "CaM", [("A1.1", 100.0, "Aspergillus flavus strain NRRL 1957 calmodulin gene"),
                                      ("A2.1", 99.9, "Aspergillus oryzae strain NRRL 447 calmodulin gene"),
                                      ("A3.1", 98.0, "Aspergillus parasiticus strain CBS 100926 calmodulin gene")])
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(flavus, 99.5)]), loci={"CaM": cam})
    assert sp == ["Aspergillus flavus", "medium"] and doc["method"] == "ITS+CaM(tie)"
    assert doc["loci"]["CaM"]["level"] == "tie" and doc["loci"]["CaM"]["candidates"] == ["Aspergillus flavus", "Aspergillus oryzae"]
    assert "tie:CaM=Aspergillus flavus/Aspergillus oryzae" in doc["flags"] and "secondary_ties_include_its_species" in doc["flags"]
    # a second, unambiguous locus resolves it
    tef = locus_b6(tmp_path, "TEF1", [("T1.1", 99.7, "Aspergillus flavus strain NRRL 1957 translation elongation factor 1-alpha gene"),
                                       ("T2.1", 99.0, "Aspergillus oryzae strain NRRL 447 translation elongation factor 1-alpha gene")])
    doc2, sp2 = run(tmp_path, helpers, unite=b6(tmp_path, [(flavus, 99.5)]), loci={"CaM": cam, "TEF1": tef})
    assert sp2 == ["Aspergillus flavus", "high"] and doc2["method"] == "ITS+TEF1(concordant)"


def test_secondary_only_when_its_is_missing(tmp_path, helpers):
    cam = locus_b6(tmp_path, "CaM", [("MH1.1", 99.8, "Aspergillus fumigatus strain CBS 133.61 calmodulin (CaM) gene")])
    tef = locus_b6(tmp_path, "TEF1", [("T1.1", 99.6, "Aspergillus fumigatus strain CBS 133.61 elongation factor 1-alpha gene")])
    doc, sp = run(tmp_path, helpers, its_text="", loci={"CaM": cam, "TEF1": tef})
    assert sp == ["Aspergillus fumigatus", "medium"] and doc["method"] == "CaM+TEF1(secondary-only)" and "no_its_species" in doc["flags"]
    doc1, sp1 = run(tmp_path, helpers, its_text="", loci={"CaM": cam})
    assert sp1 == ["Aspergillus fumigatus", "low"]


def test_medium_confidence_between_94_and_98_5(tmp_path, helpers):
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 96.0)]))
    assert sp == ["Aspergillus fumigatus", "medium"] and "below-species-threshold" in doc["method"]
    cam = locus_b6(tmp_path, "CaM", [("MH1.1", 99.8, "Aspergillus fumigatus strain CBS 133.61 calmodulin (CaM) gene")])
    doc2, sp2 = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 96.0)]), loci={"CaM": cam})
    assert sp2 == ["Aspergillus fumigatus", "medium"] and doc2["method"] == "ITS+CaM(concordant)"   # weak ITS + one locus: medium


def test_genus_only_below_94_or_sp_hit(tmp_path, helpers):
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 90.0)]))
    assert sp == ["Aspergillus sp.", "low"] and "genus-only" in doc["method"]
    doc2, sp2 = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE_SP, 99.0)]))
    assert sp2 == ["Aspergillus sp.", "low"]


def test_unknown_without_hits_and_empty_its(tmp_path, helpers):
    doc, sp = run(tmp_path, helpers, its_text="")
    assert sp == ["unknown", "none"] and doc["method"] == "no_reference" and doc["its_present"] is False and doc["its_identity"] is None


def test_gather_counts_as_a_secondary_line(tmp_path, helpers):
    g = tmp_path / "g.csv"
    g.write_text("intersect_bp,name\n1000,\"GCF_1 Aspergillus fumigatus Af293\"\n")
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 99.0)]), gather=str(g))
    assert sp == ["Aspergillus fumigatus", "high"] and doc["method"] == "ITS+genome(concordant)"
    doc2, sp2 = run(tmp_path, helpers, its_text="", gather=str(g))
    assert sp2 == ["Aspergillus fumigatus", "low"] and doc2["method"] == "genome(secondary-only)"
    g.write_text("intersect_bp,name\n1000,\"GCF_2 Candida albicans SC5314\"\n")
    doc3, sp3 = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 99.0)]), gather=str(g))
    assert sp3 == ["Aspergillus sp.", "low"] and doc3["flags"] == ["discordant:genome=Candida albicans"]


def test_species_from_title_handles_sp_brackets_and_junk():
    assert idc.species_from_title("Aspergillus flavus strain CBS 100927 calmodulin (CaM) gene, partial cds") == ("Aspergillus flavus", "Aspergillus")
    assert idc.species_from_title("Aspergillus sp. CBS 1 calmodulin gene") == (None, "Aspergillus")
    assert idc.species_from_title("[Candida] auris strain CBS 10913 28S ribosomal RNA gene") == ("Candida auris", "Candida")
    assert idc.species_from_title("uncultured fungus clone") == (None, None)
    assert idc.species_from_title("UNVERIFIED: [Neocosmospora] lechatii strain CPC 42648 cmdA gene") == ("Neocosmospora lechatii", "Neocosmospora")
    assert idc.species_from_title("TPA_asm: Fusarium solani chromosome 1") == ("Fusarium solani", "Fusarium")
    # BLAST stitle starts with the accession
    assert idc.species_from_title("AY689334.1 Aspergillus fumigatus strain KACC 41143 calmodulin gene, partial cds") == ("Aspergillus fumigatus", "Aspergillus")
    assert idc.species_from_title("NG_069679.1 Aspergillus neoglaber CBS 111.55 28S rRNA gene") == ("Aspergillus neoglaber", "Aspergillus")


def test_lsu_is_genus_level_evidence_only(tmp_path, helpers):
    lsu = locus_b6(tmp_path, "LSU", [("OL1.1", 99.9, "OL711746.1 Aspergillus spinosus strain CBS 483.65 large subunit ribosomal RNA gene", 900, 908)])
    c = idc.parse_locus_b6(lsu, "LSU")
    assert c["level"] == "genus" and c["species"] is None and c["candidates"] == ["Aspergillus spinosus"]
    doc, sp = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 99.0)]), loci={"LSU": lsu})
    assert sp == ["Aspergillus fumigatus", "medium"] and "no_secondary_species_call" in doc["flags"] and doc["flags"] == ["no_secondary_species_call"]
    assert idc.species_from_title("") == (None, None)


def test_parse_locus_b6_levels_and_coverage_filter(tmp_path):
    p = locus_b6(tmp_path, "BenA", [("A.1", 96.0, "Aspergillus terreus strain X beta-tubulin gene"),
                                     ("B.1", 99.5, "Aspergillus flavus strain Y beta-tubulin gene", 100, 520)])   # short alignment, scov 0.19
    c = idc.parse_locus_b6(p, "BenA")
    assert c["level"] == "genus" and c["genus"] == "Aspergillus" and c["species"] is None and c["n_hits"] == 1   # low-coverage hit dropped
    assert idc.parse_locus_b6(str(tmp_path / "missing.b6"), "BenA") is None
    (tmp_path / "empty.b6").write_text("")
    assert idc.parse_locus_b6(str(tmp_path / "empty.b6"), "BenA") is None
    p2 = locus_b6(tmp_path, "RPB2", [("C.1", 90.0, "Fusarium oxysporum strain Z RNA polymerase II second largest subunit gene")])
    assert idc.parse_locus_b6(p2, "RPB2")["level"] == "none"


def test_parse_mlst_and_scheme_genus_check(tmp_path, helpers):
    p = tmp_path / "mlst.tsv"
    p.write_text("asm.fa\tafumigatus\t26\tANX4(3)\tBGT1(1)\tCAT1(2)\tLIP(3)\tMAT1_2(1)\tSODB(2)\tZRF2(1)\n")
    m = idc.parse_mlst(str(p))
    assert m["scheme"] == "afumigatus" and m["st"] == "26" and m["alleles"]["MAT1_2"] == "1" and m["n_loci"] == 7 and m["n_alleles_found"] == 7
    p.write_text("asm.fa\t-\t-\n")
    assert idc.parse_mlst(str(p)) == {"scheme": None, "st": None, "alleles": {}, "n_loci": 0, "n_alleles_found": 0}
    assert idc.parse_mlst(str(tmp_path / "none.tsv")) is None
    p.write_text("asm.fa\tcalbicans\t-\tAAT1a(2)\tACC1(-)\n")
    doc, _ = run(tmp_path, helpers, unite=b6(tmp_path, [(UNITE, 99.0)]), mlst=str(p))
    assert doc["mlst"]["scheme"] == "calbicans" and doc["mlst"]["n_alleles_found"] == 1
    assert "mlst_scheme_genus_mismatch:calbicans" in doc["flags"]


def test_choose_lineage_species_then_genus_then_default(tmp_path):
    mp = tmp_path / "map.tsv"
    mp.write_text("# comment\ntaxon\tlineage\trank\nAspergillus\teurotiales_odb10\tgenus\nAspergillus niger\tspecial_odb10\tspecies\n")
    assert idc.choose_lineage("Aspergillus niger", "Aspergillus", str(mp)) == ("special_odb10", "species:Aspergillus niger")
    assert idc.choose_lineage("Aspergillus flavus", "Aspergillus", str(mp)) == ("eurotiales_odb10", "genus:Aspergillus")
    assert idc.choose_lineage("Fusarium solani", "Fusarium", str(mp)) == ("fungi_odb10", "default")
    assert idc.choose_lineage(None, None, None) == ("fungi_odb10", "no map")


def test_bundled_lineage_map_is_well_formed():
    rows = [l.rstrip("\n").split("\t") for l in open(LINEAGE_MAP) if not l.startswith("#") and l.strip()]
    assert rows[0] == ["taxon", "lineage", "rank"]
    assert all(len(r) == 3 and r[1].endswith("_odb10") and r[2] in ("genus", "species") for r in rows[1:])
    assert len({r[0] for r in rows[1:]}) == len(rows) - 1
    assert idc.choose_lineage("Candida auris", "Candida", LINEAGE_MAP)[0] == "saccharomycetes_odb10"


# ---------------------------------------------------------------- DF-005 regression (2026-09-25)
# The real CaM and RPB2 BLAST tables of ASSARM-PHI-DF-005 on the dev deployment (results_v020):
# (accession, pident, title, alen, slen). Bitscore in the fixture scales with alen, as in the run.
DF005_CAM = [
    ("LC589316.1", 99.464, "Aspergillus flavus NRRL 1957 CaM gene for calmodulin, partial cds", 746, 745),
    ("EF661506.1", 99.456, "Aspergillus oryzae isolate NRRL 447 calmodulin gene, partial cds", 735, 734),
    ("EF661508.1", 99.456, "Aspergillus flavus isolate NRRL 1957 calmodulin gene, partial cds", 735, 734),
    ("EF661516.1", 97.693, "Aspergillus parasiticus isolate NRRL 502 calmodulin gene, partial cds", 737, 736),
    ("KJ175550.1", 97.531, "Aspergillus sojae isolate CBS 100928 calmodulin (cmdA) gene, partial cds", 729, 737),
    ("EF202069.1", 99.811, "Aspergillus kambarensis strain CBS 542.69 calmodulin gene, partial cds", 530, 531),
    ("MN987053.1", 99.668, "Aspergillus agricola isolate NRRL 66869 calmodulin (cmdA) gene, partial cds", 602, 1861),
    ("EF202070.1", 99.435, "Aspergillus thomii strain CBS 120.51 calmodulin gene, partial cds", 531, 531),
]
DF005_RPB2 = [
    ("XM_032043900.1", 89.103, "Aspergillus alliaceus DNA-directed RNA polymerase II core subunit RPB2 (BDW43DRAFT_294059), mRNA", 3845, 4449),
    ("MT211766.1", 86.872, "Aspergillus hancockii strain FRR 3425 RPB2 (RPB2) gene, complete cds", 3839, 3829),
    ("MG517799.1", 99.396, "Aspergillus minisclerotigenes strain DTO 009-F7 RNA polymerase II second largest subunit (RPB2) gene, partial cds", 993, 993),
    ("MG517884.1", 99.392, "Aspergillus cerealis strain DTO 228-E7 RNA polymerase II second largest subunit (RPB2) gene, partial cds", 987, 987),
    ("MG517893.1", 99.296, "Aspergillus austwickii strain DTO 228-F7 RNA polymerase II second largest subunit (RPB2) gene, partial cds", 995, 995),
]
SYNONYMS = "name\taccepted\tbasis\nAspergillus kambarensis\tAspergillus flavus\tFrisvad 2019\nAspergillus thomii\tAspergillus flavus\tFrisvad 2019\n"


def synonyms_file(tmp_path):
    p = tmp_path / "syn.tsv"; p.write_text(SYNONYMS); return str(p)


def test_df005_cam_short_high_identity_record_widens_to_a_tie_not_a_new_species(tmp_path):
    """The 531-bp A. kambarensis record (99.81 %) beat the 746-bp A. flavus neotype (99.46 %) on identity
    alone and became the ONLY candidate, so the locus contradicted ITS. The anchor is the longest
    species-level alignment; the short record can only join the candidates."""
    call = idc.parse_locus_b6(locus_b6(tmp_path, "CaM", DF005_CAM), "CaM")
    assert call["best_accession"] == "LC589316.1" and call["best_species"] == "Aspergillus flavus"
    assert call["pident"] == 99.46 and call["top_pident"] == 99.81
    assert call["level"] == "tie" and "Aspergillus flavus" in call["candidates"] and "Aspergillus kambarensis" in call["candidates"]
    assert "Aspergillus agricola" not in call["candidates"]          # 602/1861 bp: below the subject-coverage floor


def test_df005_cam_synonyms_fold_retired_names_into_the_accepted_species(tmp_path):
    syn = idc.load_synonyms(synonyms_file(tmp_path))
    call = idc.parse_locus_b6(locus_b6(tmp_path, "CaM", DF005_CAM), "CaM", synonyms=syn)
    assert call["candidates"] == ["Aspergillus flavus", "Aspergillus oryzae"] and call["level"] == "tie"
    assert "Aspergillus kambarensis->Aspergillus flavus" in call["synonyms_applied"]
    assert "Aspergillus thomii->Aspergillus flavus" in call["synonyms_applied"]


def test_df005_rpb2_long_distant_mrna_does_not_outrank_type_material_fragments(tmp_path):
    """A 3.8-kb genome mRNA of A. alliaceus at 89 % has the top bitscore; the 1-kb section Flavi type
    records at 99.4 % carry the species signal. Anchoring on bitscore over ALL hits would call the
    locus genus-level; anchoring within the species-threshold band keeps the tie."""
    call = idc.parse_locus_b6(locus_b6(tmp_path, "RPB2", DF005_RPB2), "RPB2")
    assert call["level"] == "tie" and call["best_accession"].startswith("MG517") and call["pident"] >= 99.29
    assert set(call["candidates"]) == {"Aspergillus minisclerotigenes", "Aspergillus cerealis", "Aspergillus austwickii"}


def test_below_threshold_locus_still_reports_genus_from_the_highest_identity(tmp_path):
    rows = [("X1.1", 90.4, "Aspergillus hancockii strain FRR 3425 beta-tubulin (benA) gene, complete cds", 1597, 1926),
            ("X2.1", 96.6, "Aspergillus brasiliensis CBS 101740 beta-tubulin (ASPBRDRAFT_41588), mRNA", 900, 1758)]
    call = idc.parse_locus_b6(locus_b6(tmp_path, "BenA", rows), "BenA")
    assert call["level"] == "genus" and call["genus"] == "Aspergillus" and call["species"] is None


def test_df005_end_to_end_is_a_flavus_with_high_confidence(tmp_path, helpers):
    """ITS 100 % A. flavus, genome A. flavus, CaM a flavus/oryzae tie, RPB2 a section Flavi tie: the
    run reported `Aspergillus sp.` (low, discordant CaM=A. kambarensis). It must be A. flavus."""
    flavus = "Aspergillus_flavus|X|SH1264429.10FU|refs|k__Fungi;g__Aspergillus;s__Aspergillus_flavus"
    gather = tmp_path / "gather.csv"
    gather.write_text("intersect_bp,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund,name,filename,md5\n"
                      "30000000,0.9,0.9,0.9,0.9,1,1,0,\"GCA_002217635.1 Aspergillus flavus strain=NRRL 21882, ASM221763v1\",x,y\n")
    loci = {"CaM": locus_b6(tmp_path, "CaM", DF005_CAM), "RPB2": locus_b6(tmp_path, "RPB2", DF005_RPB2)}
    cmd_extra = ["--synonyms", synonyms_file(tmp_path)]
    its = tmp_path / "its.fa"; its.write_text(">its1\nACGT\n")
    cmd = [sys.executable, os.path.join(helpers["BIN"], "id_classify.py"), "--sample", "DF005", "--its", str(its),
           "--unite-b6", b6(tmp_path, [(flavus, 100.0)]), "--gather", str(gather),
           "--out-species", str(tmp_path / "sp.txt"), "--out-json", str(tmp_path / "id.json"), *cmd_extra]
    for name, path in loci.items():
        cmd += ["--locus-b6", f"{name}={path}"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    doc = json.load(open(tmp_path / "id.json"))
    assert doc["species"] == "Aspergillus flavus" and doc["confidence"] == "high"
    assert doc["concordance"]["secondary_disagree"] == [] and "CaM" in doc["concordance"]["ties_including_its"]
    assert not any(f.startswith("discordant") for f in doc["flags"])
    assert "tie:CaM=Aspergillus flavus/Aspergillus oryzae" in doc["flags"]


def test_load_synonyms_tolerates_comments_header_and_absence(tmp_path):
    p = tmp_path / "s.tsv"; p.write_text("# c\nname\taccepted\tbasis\nA b\tA c\tref\n\n")
    assert idc.load_synonyms(str(p)) == {"a b": "A c"} and idc.load_synonyms(None) == {} and idc.load_synonyms(str(tmp_path / "no")) == {}
    assert idc.accepted_name("A b", {"a b": "A c"}) == "A c" and idc.accepted_name("Z z", {}) == "Z z" and idc.accepted_name(None, {}) is None
