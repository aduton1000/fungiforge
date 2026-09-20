"""Unit tests for bin/validate_samplesheet.py (W4.1) — the pre-flight samplesheet check."""
import json, os, subprocess, sys
import validate_samplesheet as vs  # noqa: E402  (sys.path set in conftest)

HEADER = "sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season\n"


def sheet(tmp_path, rows, header=HEADER, name="s.csv"):
    p = tmp_path / name
    p.write_text(header + "".join(rows))
    return str(p)


def reads(tmp_path, *names):
    out = []
    for n in names:
        f = tmp_path / n
        f.write_bytes(b"@r1\nACGT\n+\nIIII\n")
        out.append(n)
    return out


def test_valid_sheet_of_three_modes(tmp_path):
    reads(tmp_path, "a.ont.fastq.gz", "b.ont.fastq.gz", "b_R1.fastq.gz", "b_R2.fastq.gz", "c_R1.fastq.gz", "c_R2.fastq.gz")
    s = sheet(tmp_path, ["A,a.ont.fastq.gz,,,AIR,F1,Dry\n",
                         "B,b.ont.fastq.gz,b_R1.fastq.gz,b_R2.fastq.gz,AIR,F1,Dry\n",
                         "C,,c_R1.fastq.gz,c_R2.fastq.gz,HUMAN,F1,Dry\n"])
    rep = vs.validate(s)
    assert rep["valid"] and rep["errors"] == [] and rep["n_samples"] == 3
    assert rep["modes"] == {"longread": 1, "hybrid": 1, "shortread": 1}
    assert rep["metadata_values"]["compartment"] == ["AIR", "HUMAN"] and rep["warnings"] == []


def test_every_error_is_reported_not_just_the_first(tmp_path):
    reads(tmp_path, "a.ont.fastq.gz", "d_R1.fastq.gz")
    (tmp_path / "empty.fastq.gz").write_bytes(b"")
    s = sheet(tmp_path, ["A,a.ont.fastq.gz,,,AIR,F1,Dry\n",
                         "A,a.ont.fastq.gz,,,AIR,F1,Dry\n",              # duplicate
                         "bad id,a.ont.fastq.gz,,,AIR,F1,Dry\n",          # unsafe id
                         ",a.ont.fastq.gz,,,AIR,F1,Dry\n",                # empty id
                         "D,,d_R1.fastq.gz,,AIR,F1,Dry\n",                # unpaired Illumina
                         "E,,,,AIR,F1,Dry\n",                             # no reads at all
                         "F,missing.fastq.gz,,,AIR,F1,Dry\n",             # missing file
                         "G,empty.fastq.gz,,,AIR,F1,Dry\n",               # empty file
                         "H,,d_R1.fastq.gz,d_R1.fastq.gz,AIR,F1,Dry\n"])  # same file twice
    rep = vs.validate(s)
    j = " | ".join(rep["errors"])
    assert not rep["valid"]
    for probe in ("duplicate sample id 'A'", "sample id 'bad id' has characters", "empty sample id",
                  "Illumina must be paired", "needs ont_fastq", "not found: 'missing.fastq.gz'",
                  "is empty: 'empty.fastq.gz'", "illumina_r1 and illumina_r2 are the same file"):
        assert probe in j, probe


def test_header_problems_and_unknown_columns(tmp_path):
    s = sheet(tmp_path, ["x,y\n"], header="id,reads\n")
    rep = vs.validate(s)
    assert not rep["valid"] and "required column 'sample' is missing" in rep["errors"][0]
    assert any("none of the read columns" in e for e in rep["errors"])
    reads(tmp_path, "a.ont.fastq.gz")
    s2 = sheet(tmp_path, ["A,a.ont.fastq.gz,,,AIR,F1,Dry,extra\n"], header=HEADER.rstrip("\n") + ",notes\n", name="s2.csv")
    rep2 = vs.validate(s2)
    assert rep2["valid"] and "columns fungiforge ignores: notes" in rep2["warnings"][0]


def test_metadata_typos_are_warned(tmp_path):
    reads(tmp_path, "a.ont.fastq.gz", "b.ont.fastq.gz")
    s = sheet(tmp_path, ["A,a.ont.fastq.gz,,,AIR,F1,Dry\n", "B,b.ont.fastq.gz,,,Air,F1 ,Dry\n"])
    rep = vs.validate(s)
    assert rep["valid"]
    j = " | ".join(rep["warnings"])
    assert "compartment: values differing only in case or spacing: AIR/Air" in j or "Air/AIR" in j
    s2 = sheet(tmp_path, ["A,a.ont.fastq.gz,,,,,\n"], name="s3.csv")
    rep2 = vs.validate(s2)
    assert rep2["valid"] and sum("empty for every isolate" in w for w in rep2["warnings"]) == 3


def test_paths_resolve_next_to_the_samplesheet(tmp_path):
    sub = tmp_path / "data"; sub.mkdir()
    (sub / "a.ont.fastq.gz").write_bytes(b"@r\nACGT\n+\nIIII\n")
    p = sub / "sheet.csv"
    p.write_text(HEADER + "A,a.ont.fastq.gz,,,AIR,F1,Dry\n")
    assert vs.validate(str(p))["valid"]                      # relative to the sheet, not the cwd
    assert vs.resolve("/abs/x.fq", str(sub)) == "/abs/x.fq"


def test_cli_and_strict_mode(tmp_path, helpers):
    reads(tmp_path, "a.ont.fastq.gz")
    s = sheet(tmp_path, ["A,a.ont.fastq.gz,,,AIR,F1,\n"])
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "validate_samplesheet.py"), "--samplesheet", s, "--json", str(tmp_path / "v.json")], capture_output=True, text=True)
    assert r.returncode == 0 and "1 isolates (1 longread): valid" in r.stdout
    assert json.load(open(tmp_path / "v.json"))["modes"] == {"longread": 1}
    r2 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "validate_samplesheet.py"), "--samplesheet", s, "--strict"], capture_output=True, text=True)
    assert r2.returncode == 1                                 # the empty season is a warning; --strict fails
    r3 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "validate_samplesheet.py"), "--samplesheet", str(tmp_path / "nope.csv")], capture_output=True, text=True)
    assert r3.returncode == 1 and "not found" in r3.stderr
    # --no-check-files validates structure without touching the reads
    s4 = sheet(tmp_path, ["A,/nowhere/a.fastq.gz,,,AIR,F1,Dry\n"], name="s4.csv")
    r4 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "validate_samplesheet.py"), "--samplesheet", s4, "--no-check-files"], capture_output=True, text=True)
    assert r4.returncode == 0


def test_bundled_test_samplesheet_is_valid(helpers):
    rep = vs.validate(os.path.join(helpers["ROOT"], "test", "samplesheet.test.csv"))
    assert rep["valid"] and rep["modes"] == {"longread": 1, "hybrid": 1, "shortread": 1}


# ── optional RNA-seq evidence columns (W6.2) ─────────────────────────────────
# rna_r1 / rna_r2 feed gene prediction. They are optional and paired: a sheet without them must
# stay valid, so every samplesheet written before this change keeps working untouched.

RNA_HEADER = "sample,ont_fastq,illumina_r1,illumina_r2,rna_r1,rna_r2,compartment,facility,season\n"


def test_a_sheet_without_rna_columns_is_unchanged(tmp_path):
    reads(tmp_path, "a.ont.fastq.gz")
    rep = vs.validate(sheet(tmp_path, ["A,a.ont.fastq.gz,,,AIR,F1,Dry\n"]))
    assert rep["valid"] and rep["n_with_rna"] == 0


def test_paired_rna_is_accepted_and_counted(tmp_path):
    reads(tmp_path, "a.ont.fastq.gz", "a_rna_R1.fastq.gz", "a_rna_R2.fastq.gz")
    s = sheet(tmp_path, ["A,a.ont.fastq.gz,,,a_rna_R1.fastq.gz,a_rna_R2.fastq.gz,AIR,F1,Dry\n"],
              header=RNA_HEADER)
    rep = vs.validate(s)
    assert rep["valid"], rep["errors"]
    assert rep["n_with_rna"] == 1


def test_unpaired_rna_is_refused(tmp_path):
    reads(tmp_path, "a.ont.fastq.gz", "a_rna_R1.fastq.gz")
    s = sheet(tmp_path, ["A,a.ont.fastq.gz,,,a_rna_R1.fastq.gz,,AIR,F1,Dry\n"], header=RNA_HEADER)
    rep = vs.validate(s)
    assert not rep["valid"]
    assert any("RNA-seq must be paired" in e for e in rep["errors"])


def test_the_same_file_twice_is_refused(tmp_path):
    reads(tmp_path, "a.ont.fastq.gz", "a_rna_R1.fastq.gz")
    s = sheet(tmp_path, ["A,a.ont.fastq.gz,,,a_rna_R1.fastq.gz,a_rna_R1.fastq.gz,AIR,F1,Dry\n"],
              header=RNA_HEADER)
    rep = vs.validate(s)
    assert not rep["valid"]
    assert any("same file" in e for e in rep["errors"])


def test_a_missing_rna_file_is_caught_at_the_sheet(tmp_path):
    """The point of the pre-flight: fail here, not four hours later inside prediction."""
    reads(tmp_path, "a.ont.fastq.gz", "a_rna_R1.fastq.gz")
    s = sheet(tmp_path, ["A,a.ont.fastq.gz,,,a_rna_R1.fastq.gz,gone_R2.fastq.gz,AIR,F1,Dry\n"],
              header=RNA_HEADER)
    rep = vs.validate(s)
    assert not rep["valid"]
    assert any("rna_r2 not found" in e for e in rep["errors"])


def test_rna_columns_are_not_reported_as_unknown(tmp_path):
    reads(tmp_path, "a.ont.fastq.gz", "a_rna_R1.fastq.gz", "a_rna_R2.fastq.gz")
    s = sheet(tmp_path, ["A,a.ont.fastq.gz,,,a_rna_R1.fastq.gz,a_rna_R2.fastq.gz,AIR,F1,Dry\n"],
              header=RNA_HEADER)
    rep = vs.validate(s)
    assert not any("unknown" in w.lower() for w in rep["warnings"]), rep["warnings"]
