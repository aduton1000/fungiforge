"""Unit tests for bin/validate_master.py and bin/cohort_report.py (W3.3): the master-table schema,
the merged table, the cohort summary and the MultiQC custom content."""
import csv, json, os, subprocess, sys
import validate_master as vm  # noqa: E402  (sys.path set in conftest)
import cohort_report as cr  # noqa: E402
import make_report as mr  # noqa: E402

SCHEMA = os.path.join(os.path.dirname(__file__), "..", "..", "fungiforge", "resources", "master_schema.json")


def test_schema_matches_the_master_columns_exactly():
    s = vm.load_schema(SCHEMA)
    names = [c["name"] for c in s["columns"]]
    assert names == mr.MASTER_COLS, "master_schema.json and MASTER_COLS have drifted apart"
    assert s["primary_key"] == "sample" and s["na_value"] == "NA" and s["delimiter"] == "\t"
    assert all(c["type"] in ("string", "integer", "number", "boolean") and c.get("description") and c.get("group") for c in s["columns"])
    assert {c["name"] for c in s["columns"] if c.get("enum")} >= {"species_confidence", "novelty", "polish_mode", "mating_type"}


def write_master(path, rows, cols=None):
    cols = cols or mr.MASTER_COLS
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "NA") for c in cols})
    return str(path)


def base_row(sample="S1", **kw):
    r = {c: "NA" for c in mr.MASTER_COLS}
    r.update(sample=sample, species="Aspergillus fumigatus", species_confidence="high", qc_pass="True",
             busco_complete="99.7", assembly_len="29986658", n_contigs="48", resistant_classes="none",
             n_known_af_mutations="0", novelty="known_species", gate="pass", stages_failed="none", polish_mode="hybrid")
    r.update(kw)
    return r


def test_validate_master_accepts_a_good_table_and_reports_every_error(tmp_path):
    schema = vm.load_schema(SCHEMA)
    good = write_master(tmp_path / "good.tsv", [base_row("S1"), base_row("S2", compartment="AIR")])
    rep = vm.validate(good, schema)
    assert rep["valid"] and rep["n_rows"] == 2 and rep["errors"] == [] and rep["warnings"] == []
    bad = write_master(tmp_path / "bad.tsv", [base_row("S1", busco_complete="high", n_contigs="4.5", qc_pass="maybe",
                                                       species_confidence="excellent", polish_mode="nanopore"),
                                              base_row("S1")])
    rep2 = vm.validate(bad, schema)
    joined = " | ".join(rep2["errors"])
    assert not rep2["valid"]
    for probe in ("busco_complete: 'high' is not a number", "n_contigs: '4.5' is not an integer", "qc_pass: 'maybe' is not a boolean",
                  "species_confidence: 'excellent' is not one of", "polish_mode: 'nanopore' is not one of", "duplicate sample 'S1'"):
        assert probe in joined, probe


def test_validate_master_missing_and_extra_columns(tmp_path):
    schema = vm.load_schema(SCHEMA)
    short = write_master(tmp_path / "short.tsv", [base_row()], cols=mr.MASTER_COLS[:-3])
    rep = vm.validate(short, schema)
    assert not rep["valid"] and "missing columns" in rep["errors"][0]
    extra = write_master(tmp_path / "extra.tsv", [dict(base_row(), future_col="x")], cols=mr.MASTER_COLS + ["future_col"])
    rep2 = vm.validate(extra, schema)
    assert rep2["valid"] and "future_col" in rep2["warnings"][0]        # append-only: newer tables still validate


def test_validate_master_cli(tmp_path, helpers):
    good = write_master(tmp_path / "m.tsv", [base_row()])
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "validate_master.py"), "--master", good, "--json", str(tmp_path / "v.json")], capture_output=True, text=True)
    assert r.returncode == 0 and "valid" in r.stdout and json.load(open(tmp_path / "v.json"))["n_schema_columns"] == len(mr.MASTER_COLS)
    bad = write_master(tmp_path / "b.tsv", [base_row(n_contigs="many")])
    r2 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "validate_master.py"), "--master", bad], capture_output=True, text=True)
    assert r2.returncode == 1 and "not an integer" in r2.stderr


def test_read_masters_merges_and_orders(tmp_path):
    cols = mr.MASTER_COLS
    write_master(tmp_path / "S2.master.tsv", [base_row("S2")])
    write_master(tmp_path / "S1.master.tsv", [base_row("S1")])
    rows = cr.read_masters([str(tmp_path / "S2.master.tsv"), str(tmp_path / "S1.master.tsv")], cols)
    assert [r["sample"] for r in rows] == ["S1", "S2"] and list(rows[0]) == cols


def test_summarise_counts_and_metrics():
    rows = [base_row("S1", compartment="AIR", resistant_classes="azole", mycotoxin_clusters="gliotoxin", n_bgc="52", busco_complete="99.7"),
            base_row("S2", compartment="AIR", resistant_classes="azole;echinocandin", mycotoxin_clusters="aflatoxin;sterigmatocystin", n_bgc="102", busco_complete="99.1", species="Aspergillus flavus"),
            base_row("S3", compartment="HUMAN", gate="skipped(read_triage:non_fungal)", species="unknown", stages_failed="gate:skipped", n_bgc="NA", busco_complete="NA")]
    s = cr.summarise(rows)
    assert s["n_isolates"] == 3 and s["n_with_resistance"] == 2 and s["n_gated"] == 1 and s["n_stage_failures"] == 1
    assert dict(s["species"])["Aspergillus fumigatus"] == 1 and dict(s["compartments"])["AIR"] == 2
    assert dict(s["resistant_classes"]) == {"azole": 2, "echinocandin": 1}
    assert dict(s["mycotoxins"]) == {"gliotoxin": 1, "aflatoxin": 1, "sterigmatocystin": 1}
    assert s["metrics"]["n_bgc"] == {"n": 2, "min": 52.0, "median": 77.0, "max": 102.0, "mean": 77.0}
    assert s["metrics"]["busco_complete"]["n"] == 2
    s2 = cr.summarise(rows, cohort={"n_isolates": 3, "species_clusters": [{"cluster": "C1", "members": ["S1"], "species": ["Aspergillus fumigatus"]}], "phylogenomics": {"tree": None}},
                      provenance={"pipeline": {"version": "0.2.0", "commit_id": "abc"}, "run": {"run_name": "r1", "start": "2026"}})
    assert s2["cohort"]["species_clusters"][0]["cluster"] == "C1" and s2["provenance"]["pipeline"] == "0.2.0"


def test_cohort_report_cli_writes_merged_table_report_and_multiqc(tmp_path, helpers):
    md = tmp_path / "masters"; md.mkdir()
    write_master(md / "S1.master.tsv", [base_row("S1", mycotoxin_clusters="gliotoxin")])
    write_master(md / "S2.master.tsv", [base_row("S2", species="Aspergillus flavus", resistant_classes="azole")])
    out = tmp_path / "04_summary"
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "cohort_report.py"), "--masters", str(md), "--outdir", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "2 isolates" in r.stdout and "schema valid" in r.stdout
    merged = list(csv.DictReader(open(out / "master_fungi.tsv"), delimiter="\t"))
    assert [m["sample"] for m in merged] == ["S1", "S2"] and list(merged[0]) == mr.MASTER_COLS
    s = json.load(open(out / "cohort_summary.json"))
    assert s["n_isolates"] == 2 and s["schema_validation"]["valid"] and s["n_with_resistance"] == 1
    html = (out / "cohort_report.html").read_text()
    assert "FungiForge cohort" in html and "Aspergillus flavus" in html and "gliotoxin" in html
    mq = json.load(open(out / "multiqc" / "fungiforge_mqc.json"))
    assert mq["plot_type"] == "table" and set(mq["data"]) == {"S1", "S2"} and mq["data"]["S1"]["busco_complete"] == 99.7
    assert (out / "multiqc" / "multiqc_config.yaml").exists() and (out / "multiqc" / "fungiforge_species_mqc.json").exists()
    # an invalid row makes the stage fail loudly, and the merged table is still written for inspection
    write_master(md / "S3.master.tsv", [base_row("S3", n_contigs="lots")])
    r2 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "cohort_report.py"), "--masters", str(md), "--outdir", str(out), "--no-multiqc-export"], capture_output=True, text=True)
    assert r2.returncode == 1 and "INVALID" in r2.stdout and (out / "master_fungi.tsv").exists()
