"""Unit tests for bin/make_report.py — master-table row construction and the status columns."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "bin")))
import make_report as m  # noqa: E402


def row(S):
    return m.build_row("X", "AIR", "F1", "Dry", S)


def test_columns_are_append_only_and_end_with_status_columns():
    assert m.MASTER_COLS[:4] == ["sample", "compartment", "facility", "season"]
    assert m.MASTER_COLS[-5:] == ["mito_size_kb", "mito_core_genes", "mito_circular", "mito_copy_ratio", "mito_heteroplasmic_sites"]
    assert m.MASTER_COLS[-9:-5] == ["n_secreted", "n_effectors", "n_cazymes", "n_phibase_hits"]
    assert m.MASTER_COLS[-32:-9] == ["sample_verdict", "contam_removed_pct", "top_taxon", "stages_failed", "gate",
                                   "read_verdict", "genome_size_est", "heterozygosity_pct", "ploidy_hint", "coverage",
                                   "id_loci_agree", "id_flags", "mlst_st", "busco_lineage_specific", "busco_complete_specific",
                                   "resistance_read_support", "copy_number_flags",
                                   "n_proteins", "pct_pfam", "pct_go", "pct_eggnog", "pct_interpro", "annotation_training"]
    assert len(m.MASTER_COLS) == len(set(m.MASTER_COLS))


def test_row_has_every_column():
    r = row({})
    assert set(r) == set(m.MASTER_COLS)


def test_stages_failed_none_when_all_ok():
    S = {"annotate": {"status": "ok"}, "bgc": {"status": "ok", "n_clusters": 3}}
    assert row(S)["stages_failed"] == "none"


def test_stages_failed_lists_non_ok_stages_sorted():
    S = {"bgc": {"status": "failed"}, "polish": {"status": "partial"}, "novelty": {"status": "skipped"}, "annotate": {"status": "ok"}}
    assert row(S)["stages_failed"] == "bgc:failed;novelty:skipped;polish:partial"


def test_stages_without_status_are_not_flagged():
    S = {"readqc": {"illumina": True}}
    assert row(S)["stages_failed"] == "none"


def test_bgc_failed_gives_na_not_zero():
    assert row({"bgc": {"status": "failed", "n_clusters": None}})["n_bgc"] == "NA"
    assert row({"bgc": {"status": "ok", "n_clusters": 51}})["n_bgc"] == 51


def test_decontam_columns():
    S = {"decontam": {"verdict": "fungal", "dropped_bp_pct": 0.04,
                      "top_species": [{"name": "Aspergillus flavus", "pct": 51.3}]}}
    r = row(S)
    assert r["sample_verdict"] == "fungal" and r["contam_removed_pct"] == 0.04
    assert r["top_taxon"] == "Aspergillus flavus (51.3%)"


# ---------------------------------------------------------------- CLI end to end (W0.3)
import json, subprocess, sys  # noqa: E402


def test_cli_merges_stage_jsons_into_html_and_master_row(tmp_path):
    stages = {
        "identify": {"stage": "identify", "species": "Aspergillus fumigatus", "confidence": "high", "method": "ITS/UNITE", "status": "ok"},
        "assembly_qc": {"stage": "assembly_qc", "qc_pass": True, "busco_complete": 99.1, "lineage": "fungi_odb10", "assembly_len": 29000000, "n_contigs": 43, "status": "ok"},
        "resistance": {"stage": "resistance", "status": "ok", "polish_mode": "hybrid",
                       "calls": [{"gene": "cyp51A", "drug_class": "azole", "change": "L98H", "known": True, "confidence": "high"}],
                       "cyp51A_TR": {"tr_type": "TR34"}, "summary": {"resistant_drug_classes": ["azole"], "n_known_mutations": 1}},
        "bgc": {"stage": "bgc", "status": "failed", "n_clusters": None},
        "decontam": {"stage": "decontam", "status": "ok", "verdict": "fungal", "dropped_bp_pct": 0.5, "top_species": [{"name": "Aspergillus fumigatus", "pct": 30.0}]},
    }
    for k, v in stages.items():
        (tmp_path / f"S1.{k}.json").write_text(json.dumps(v))
    (tmp_path / "S1.broken.json").write_text("{not json")
    r = subprocess.run([sys.executable, os.path.join(HERE, "..", "..", "bin", "make_report.py"), "--sample", "S1", "--compartment", "AIR",
                        "--facility", "F1", "--season", "Dry", "--jsons", str(tmp_path / "S1.*.json"),
                        "--html", str(tmp_path / "S1.html"), "--master", str(tmp_path / "S1.master.tsv")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    header, row = [l.rstrip("\n").split("\t") for l in open(tmp_path / "S1.master.tsv")]
    assert header == m.MASTER_COLS
    d = dict(zip(header, row))
    assert d["sample"] == "S1" and d["species"] == "Aspergillus fumigatus" and d["qc_pass"] == "True"
    assert d["resistant_classes"] == "azole" and d["cyp51A_TR"] == "TR34" and d["n_bgc"] == "NA" and d["stages_failed"] == "bgc:failed"
    assert d["top_taxon"] == "Aspergillus fumigatus (30.0%)" and d["sample_verdict"] == "fungal" and d["contam_removed_pct"] == "0.5"
    html_text = open(tmp_path / "S1.html").read()
    assert "<title>FungiForge · S1</title>" in html_text and "L98H" in html_text and "Aspergillus fumigatus" in html_text
    assert "5 stages" in r.stdout          # the unreadable JSON was skipped, not fatal



def test_gate_column_reports_skipped_reason_and_stages_failed():
    S = {"decontam": {"status": "ok", "verdict": "non_fungal"},
         "gate": {"stage": "gate", "status": "skipped", "reason": "verdict:non_fungal"}}
    r = row(S)
    assert r["gate"] == "skipped(verdict:non_fungal)" and r["stages_failed"] == "gate:skipped"
    assert r["species"] == "NA" and r["n_bgc"] == "NA"
    assert row({"decontam": {"status": "ok", "verdict": "fungal"}})["gate"] == "pass"


def test_read_triage_and_kmer_columns():
    S = {"triage": {"status": "ok", "verdict": "likely_fungal"},
         "kmer": {"status": "partial", "genome_size_est": 28864332, "heterozygosity_pct": None, "ploidy_hint": "haploid",
                  "coverage_from_kmers": 40.5, "coverage_from_read_bases": 50.0}}
    r = row(S)
    assert r["read_verdict"] == "likely_fungal" and r["genome_size_est"] == 28864332 and r["heterozygosity_pct"] == "NA"
    assert r["ploidy_hint"] == "haploid" and r["coverage"] == 40.5 and r["stages_failed"] == "kmer:partial"
    r2 = row({"kmer": {"status": "ok", "genome_size_est": None, "coverage_from_kmers": None, "coverage_from_read_bases": 50.0}})
    assert r2["coverage"] == 50.0 and r2["genome_size_est"] == "NA" and r2["read_verdict"] == "NA"
    S3 = {"triage": {"status": "ok", "verdict": "non_fungal"}, "gate": {"stage": "gate", "status": "skipped", "reason": "read_triage:non_fungal"}}
    assert row(S3)["gate"] == "skipped(read_triage:non_fungal)" and row(S3)["read_verdict"] == "non_fungal"


def test_top_taxon_falls_back_to_the_read_triage_when_decontam_never_ran():
    S = {"triage": {"status": "ok", "verdict": "non_fungal", "top_species": [{"name": "Klebsiella pneumoniae", "pct": 71.2}]},
         "gate": {"stage": "gate", "status": "skipped", "reason": "read_triage:non_fungal"}}
    r = row(S)
    assert r["top_taxon"] == "Klebsiella pneumoniae (71.2%)" and r["sample_verdict"] == "NA" and r["read_verdict"] == "non_fungal"


def test_identification_columns_w23():
    S = {"identify": {"species": "Aspergillus fumigatus", "confidence": "high", "method": "ITS+CaM(concordant)",
                      "flags": ["tie:BenA=Aspergillus fumigatus/Aspergillus lentulus"],
                      "concordance": {"secondary_agree": ["CaM", "TEF1"]}, "mlst": {"scheme": "afumigatus", "st": "26"}},
         "busco_lineage": {"lineage": "eurotiales_odb10", "busco_complete": 99.4}}
    r = row(S)
    assert r["id_loci_agree"] == "CaM;TEF1" and r["id_flags"] == "tie:BenA=Aspergillus fumigatus/Aspergillus lentulus"
    assert r["mlst_st"] == "afumigatus:ST26" and r["busco_lineage_specific"] == "eurotiales_odb10" and r["busco_complete_specific"] == 99.4
    r2 = row({"identify": {"species": "x", "mlst": {"scheme": "calbicans", "st": None}}, "busco_lineage": {"lineage": "fungi_odb10", "busco_complete": None, "skipped": "same_as_qc"}})
    assert r2["mlst_st"] == "calbicans:ST-" and r2["busco_lineage_specific"] == "NA" and r2["busco_complete_specific"] == "NA"
    r3 = row({})
    assert r3["id_loci_agree"] == "none" and r3["id_flags"] == "none" and r3["mlst_st"] == "NA"


def test_resistance_read_support_columns_w24():
    S = {"resistance": {"summary": {"read_support": {"mode": "reads", "confirmed": 2, "discordant": 0, "reads_only": 1, "insufficient": 3},
                                    "copy_number_flags": ["ERG11:2.1"]}}}
    r = row(S)
    assert r["resistance_read_support"] == "confirmed:2;discordant:0;reads_only:1;insufficient:3" and r["copy_number_flags"] == "ERG11:2.1"
    assert row({"resistance": {"summary": {"read_support": {"mode": "assembly_only"}}}})["resistance_read_support"] == "assembly_only"
    assert row({})["resistance_read_support"] == "NA" and row({})["copy_number_flags"] == "none"


def test_annotation_columns_w25():
    S = {"annotate": {"n_proteins": 9607, "pct_pfam": 61.2, "pct_go": 40.0, "pct_eggnog": 70.5, "pct_interpro": 66.1},
         "predict": {"augustus_species": "aspergillus_fumigatus", "busco_db": "eurotiomycetes", "genemark": "yes"}}
    r = row(S)
    assert (r["n_proteins"], r["pct_pfam"], r["pct_eggnog"]) == (9607, 61.2, 70.5) and r["annotation_training"] == "aspergillus_fumigatus/eurotiomycetes/genemark:yes"
    assert row({})["annotation_training"] == "NA" and row({})["n_proteins"] == "NA"


def test_organelle_columns_w27():
    S = {"organelle": {"mito_present": True, "mito_size": 30696, "n_core_genes": 15, "circular": True, "copy_ratio": 35.2, "heteroplasmy": {"n_heteroplasmic_sites": 2}}}
    r = row(S)
    assert (r["mito_size_kb"], r["mito_core_genes"], r["mito_circular"], r["mito_copy_ratio"], r["mito_heteroplasmic_sites"]) == (30.7, "15/15", True, 35.2, 2)
    r2 = row({"organelle": {"mito_present": False, "mito_size": 0, "n_core_genes": 0}})
    assert (r2["mito_size_kb"], r2["mito_core_genes"], r2["mito_circular"]) == ("NA", "NA", "NA")
