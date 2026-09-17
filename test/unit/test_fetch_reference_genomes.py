"""Unit tests for bin/fetch_reference_genomes.py with the NCBI Datasets API replaced by canned responses."""
import io, json, os, zipfile
import fetch_reference_genomes as fg  # noqa: E402  (sys.path set in conftest)

GENERA = os.path.join(os.path.dirname(__file__), "..", "..", "fungiforge", "resources", "novelty_genera.txt")


def report(acc, name, strain="", cat="representative genome", length=30000000):
    return {"accession": acc, "organism": {"organism_name": name, "tax_id": 1, "infraspecific_names": {"strain": strain}},
            "assembly_info": {"refseq_category": cat, "assembly_level": "Chromosome"}, "assembly_stats": {"total_sequence_length": str(length)}}


def test_list_reference_genomes_paginates(monkeypatch):
    pages = {None: {"reports": [report("GCF_1", "Aspergillus fumigatus Af293", "Af293")], "next_page_token": "T2"},
             "T2": {"reports": [report("GCA_2", "Aspergillus lentulus", cat="")]}}
    def fake_get(url, **kw):
        tok = None
        if "page_token=" in url:
            tok = url.split("page_token=")[1].split("&")[0]
        assert "reference_only" in url and "Aspergillus" in url
        return json.dumps(pages[tok])
    monkeypatch.setattr(fg, "_get", fake_get)
    rows = fg.list_reference_genomes("Aspergillus")
    assert [r["accession"] for r in rows] == ["GCF_1", "GCA_2"] and rows[0]["species"] == "Aspergillus fumigatus" and rows[0]["strain"] == "Af293"
    assert fg.list_reference_genomes("Aspergillus", max_n=1)[0]["accession"] == "GCF_1"


def test_download_fasta_extracts_fna(tmp_path, monkeypatch):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("ncbi_dataset/data/GCF_1/GCF_1_x_genomic.fna", ">c1\nACGT\n")
        z.writestr("ncbi_dataset/data/GCF_1/README.md", "x")
    monkeypatch.setattr(fg, "_get", lambda url, **kw: buf.getvalue())
    out = tmp_path / "GCF_1.fna"
    assert fg.download_fasta("GCF_1", str(out)) == 1 and out.read_text() == ">c1\nACGT\n" and not (tmp_path / "GCF_1.fna.part").exists()


def test_accession_report_fills_metadata(monkeypatch):
    monkeypatch.setattr(fg, "_get", lambda url, **kw: json.dumps({"reports": [report("GCF_000002655.1", "Aspergillus fumigatus Af293", "Af293", "reference genome")]}))
    r = fg.accession_report("GCF_000002655.1")
    assert r["species"] == "Aspergillus fumigatus" and r["strain"] == "Af293" and r["category"] == "reference genome" and r["source"] == "accession"
    monkeypatch.setattr(fg, "_get", lambda url, **kw: (_ for _ in ()).throw(RuntimeError("down")))
    assert fg.accession_report("GCF_X")["organism"] == "" and fg.accession_report("GCF_X")["category"] == "explicit"


def test_bundled_genus_list_and_cli_dry_run(tmp_path, monkeypatch, helpers):
    genera = [l.strip() for l in open(GENERA) if l.strip() and not l.startswith("#")]
    assert "Aspergillus" in genera and "Candida" in genera and "Cryptococcus" in genera and len(genera) > 80 and genera == sorted(genera)
    import subprocess, sys
    (tmp_path / "g.txt").write_text("Aspergillus\n")
    # dry run against the network is avoided: point --genera at nothing and pass an explicit accession
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "fetch_reference_genomes.py"), "--out-dir", str(tmp_path / "db"), "--genera", ",",
                        "--accessions", "GCF_000002655.1", "--dry-run"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    man = (tmp_path / "db" / "manifest.tsv").read_text().splitlines()
    assert man[0].startswith("accession\torganism") and man[1].startswith("GCF_000002655.1\t") and "1 genomes listed" in r.stdout
