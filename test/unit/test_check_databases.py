"""Unit tests for bin/check_databases.py (DB_CHECK stage)."""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.abspath(os.path.join(HERE, "..", "..", "bin"))
SCRIPT = os.path.join(BIN, "check_databases.py")


def make_db(root, names, marker=True, manifest=True):
    """Create a fake --data_dir with the key files each named database needs."""
    files = {
        "unite": ["sh_general_release_dynamic_19.02.2025.fasta"],
        "kraken2": ["hash.k2d", "taxo.k2d", "opts.k2d"],
        "busco": ["lineages/fungi_odb10/dataset.cfg"],
        "funannotate": ["Pfam-A.hmm", "funannotate-db-info.txt"],
        "antismash": ["clusterblast/x", "knownclusterblast/x", "pfam/x"],
        "fungamr": ["FungAMR_070425.tsv", "reference_proteins.faa"],
        "refseq_fungi": ["genbank-2022.03-fungi-k31.zip"],
    }
    rows = ["database\tdetail\tstatus\ttimestamp"]
    for n in names:
        d = os.path.join(root, n)
        for f in files[n]:
            p = os.path.join(d, f); os.makedirs(os.path.dirname(p), exist_ok=True)
            open(p, "w").write("x")
        if marker:
            open(os.path.join(d, ".done"), "w").write("2026-09-01 12:00:00")
        rows.append(f"{n}\tsource-of-{n}\tOK\t2026-09-01 12:00:00")
    if manifest:
        open(os.path.join(root, "MANIFEST.tsv"), "w").write("\n".join(rows) + "\n")


def run(root, out, *extra):
    r = subprocess.run([sys.executable, SCRIPT, "--data-dir", str(root), "--out", str(out), *extra],
                       capture_output=True, text=True)
    doc = json.load(open(out)) if os.path.exists(out) else None
    return r.returncode, doc, r.stdout + r.stderr


ALL = ["unite", "kraken2", "busco", "funannotate", "antismash", "fungamr", "refseq_fungi"]


def test_complete_database_root_passes_and_records_sources(tmp_path):
    make_db(tmp_path, ALL)
    rc, doc, out = run(tmp_path, tmp_path / "m.json")
    assert rc == 0 and doc["ok"] and doc["problems"] == []
    assert doc["databases"]["unite"]["done_marker"].startswith("2026-09-01")
    assert doc["databases"]["kraken2"]["source"]["detail"] == "source-of-kraken2"
    assert doc["databases"]["busco"]["missing"] == []          # any-of layout accepted
    assert "rvdb" in doc["optional"] and doc["optional"]["rvdb"]["present"] is False


def test_missing_required_database_fails_with_named_problem(tmp_path):
    make_db(tmp_path, [n for n in ALL if n != "kraken2"])
    rc, doc, out = run(tmp_path, tmp_path / "m.json")
    assert rc == 1 and not doc["ok"]
    assert any(p.startswith("kraken2: missing directory") for p in doc["problems"])
    assert "refusing to start" in out


def test_missing_done_marker_is_reported_as_incomplete(tmp_path):
    make_db(tmp_path, ALL, marker=False)
    rc, doc, _ = run(tmp_path, tmp_path / "m.json")
    assert rc == 1 and all("no .done marker" in p for p in doc["problems"]) and len(doc["problems"]) == len(ALL)


def test_allow_missing_downgrades_to_warning(tmp_path):
    make_db(tmp_path, ["unite"])
    rc, doc, _ = run(tmp_path, tmp_path / "m.json", "--allow-missing")
    assert rc == 0 and not doc["ok"] and len(doc["problems"]) == len(ALL) - 1


def test_require_subset_only_checks_those(tmp_path):
    make_db(tmp_path, ["unite", "busco"])
    rc, doc, _ = run(tmp_path, tmp_path / "m.json", "--require", "unite,busco")
    assert rc == 0 and doc["ok"] and set(doc["databases"]) == {"unite", "busco"}


def test_path_override_honours_funannotate_db_elsewhere(tmp_path):
    make_db(tmp_path, [n for n in ALL if n != "funannotate"])
    other = tmp_path / "elsewhere"
    make_db(other, ["funannotate"], manifest=False)
    rc, doc, _ = run(tmp_path, tmp_path / "m.json", "--path", f"funannotate={other / 'funannotate'}")
    assert rc == 0 and doc["databases"]["funannotate"]["path"] == str(other / "funannotate")


def test_partial_key_files_are_listed(tmp_path):
    make_db(tmp_path, ALL)
    os.remove(tmp_path / "kraken2" / "taxo.k2d")
    rc, doc, _ = run(tmp_path, tmp_path / "m.json")
    assert rc == 1 and doc["databases"]["kraken2"]["missing"] == ["taxo.k2d"]
