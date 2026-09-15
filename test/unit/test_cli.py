"""Unit tests for the fungiforge console entry point (fungiforge/cli.py): samplesheet building
from read directories, and the run / fetch-refs / version commands."""
import os, subprocess, sys
import pytest
import fungiforge.cli as cli  # noqa: E402  (repo root on sys.path via conftest)


def touch(p):
    os.makedirs(os.path.dirname(p), exist_ok=True); open(p, "w").write("@r\nA\n+\nI\n"); return p


def sheet_rows(path):
    return [l.rstrip("\n").split(",") for l in open(path)][1:]


def test_samplesheet_routes_ont_hybrid_and_illumina_only(tmp_path, capsys):
    ont = tmp_path / "ont"; il = tmp_path / "il"
    touch(str(ont / "ISO1.ont.fastq.gz"))                                  # ONT only
    touch(str(ont / "ISO2.fastq.gz")); touch(str(il / "ISO2_R1.fastq.gz")); touch(str(il / "ISO2_R2.fastq.gz"))   # hybrid
    touch(str(il / "ISO3_S7_L001_R1_001.fastq.gz")); touch(str(il / "ISO3_S7_L001_R2_001.fastq.gz"))              # Illumina only, bcl2fastq names
    touch(str(il / "ISO4_1.fq.gz")); touch(str(il / "ISO4_2.fq.gz"))                                              # SRA-style names
    touch(str(il / "LONELY_R1.fastq.gz"))                                                                        # unpaired -> dropped
    out = tmp_path / "s.csv"
    cli.main(["samplesheet", "--ont", str(ont / "*.fastq.gz"), "--illumina-r1", str(il / "*_R1*"), str(il / "*_1.fq.gz"),
              "--illumina-r2", str(il / "*_R2*"), str(il / "*_2.fq.gz"), "--compartment", "AIR", "--facility", "F1", "--season", "Dry", "-o", str(out)])
    rows = {r[0]: r for r in sheet_rows(out)}
    assert set(rows) == {"ISO1", "ISO2", "ISO3", "ISO4"}
    assert rows["ISO1"][1].endswith("ISO1.ont.fastq.gz") and rows["ISO1"][2] == "" and rows["ISO1"][3] == ""
    assert rows["ISO2"][1].endswith("ISO2.fastq.gz") and rows["ISO2"][2].endswith("ISO2_R1.fastq.gz") and rows["ISO2"][3].endswith("ISO2_R2.fastq.gz")
    assert rows["ISO3"][1] == "" and rows["ISO3"][2].endswith("_R1_001.fastq.gz") and rows["ISO3"][3].endswith("_R2_001.fastq.gz")
    assert rows["ISO4"][2].endswith("ISO4_1.fq.gz") and rows["ISO4"][3].endswith("ISO4_2.fq.gz")
    assert rows["ISO1"][4:] == ["AIR", "F1", "Dry"]
    cap = capsys.readouterr().out
    assert "ONT-only=1, Illumina-only=2, hybrid=1" in cap and "skipped 1 sample(s)" in cap and "LONELY" in cap


def test_samplesheet_header_matches_pipeline_contract(tmp_path):
    touch(str(tmp_path / "A.ont.fastq.gz"))
    out = tmp_path / "s.csv"
    cli.main(["samplesheet", "--ont", str(tmp_path / "*.ont.fastq.gz"), "-o", str(out)])
    assert open(out).readline().strip() == "sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season"
    assert sheet_rows(out)[0][4:] == ["NA", "NA", "NA"]


def test_samplesheet_with_no_matching_reads_exits(tmp_path):
    with pytest.raises(SystemExit):
        cli.main(["samplesheet", "--ont", str(tmp_path / "nothing*.fastq.gz"), "-o", str(tmp_path / "s.csv")])


def test_run_builds_the_nextflow_command(monkeypatch, capsys):
    seen = {}
    monkeypatch.setattr(cli.subprocess, "call", lambda cmd: seen.setdefault("cmd", cmd) and 0)
    rc = cli.main(["run", "--samplesheet", "s.csv", "--data_dir", "/db", "--outdir", "res", "-profile", "hpc_slurm,apptainer", "-resume", "--", "--skip_bgc"])
    cmd = seen["cmd"]
    assert cmd[:2] == ["nextflow", "run"] and cmd[2].endswith("main.nf")
    assert cmd[3:5] == ["-profile", "hpc_slurm,apptainer"]
    assert "--samplesheet" in cmd and "--data_dir" in cmd and "--outdir" in cmd and "-resume" in cmd and cmd[-1] == "--skip_bgc"
    assert "[fungiforge] nextflow run" in capsys.readouterr().out


def test_fetch_refs_and_version(capsys):
    cli.main(["fetch-refs", "--data_dir", "/data/ff"])
    out = capsys.readouterr().out
    assert 'export FUNGIFORGE_DB="/data/ff"' in out and out.strip().endswith("fetch_references.sh")
    cli.main(["version"])
    assert capsys.readouterr().out.startswith("FungiForge v")


def test_console_script_entry_point_runs():
    r = subprocess.run([sys.executable, "-m", "fungiforge.cli", "version"], capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    assert r.returncode == 0 and r.stdout.startswith("FungiForge v")
