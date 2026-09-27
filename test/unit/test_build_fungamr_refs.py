"""Unit tests for bin/build_fungamr_refs.py — FungAMR table parsing and the UniProt fetch logic,
with the network call replaced by a fake (no network in tests)."""
import json, os, sys
import build_fungamr_refs as bf  # noqa: E402  (sys.path set in conftest)


def test_sanitize_and_fasta_parse():
    assert bf.sanitize("Aspergillus fumigatus (Af293)") == "Aspergillus_fumigatus_Af293"
    assert bf.sanitize("  cyp51A/erg11 ") == "cyp51A_erg11"
    out = {}
    bf._parse(">sp|Q4WNT5|CP51A_ASPFU Sterol 14-alpha demethylase\nMGLIA\nFVAG\n>tr|A0A0J5|X\nMKK\n", out)
    assert out == {"Q4WNT5": "MGLIAFVAG", "A0A0J5": "MKK"}


def test_fetch_fasta_splits_batches_around_bad_accessions(monkeypatch):
    calls = []
    def fake_get(batch):
        calls.append(list(batch))
        if "BAD1" in batch:
            return None                                     # UniProt 400s the whole batch
        return "".join(f">sp|{a}|N\nSEQ{a}\n" for a in batch)
    monkeypatch.setattr(bf, "_get", fake_get)
    seqs = bf.fetch_fasta(["A1", "A2", "BAD1", "A3"], chunk=80)
    assert seqs == {"A1": "SEQA1", "A2": "SEQA2", "A3": "SEQA3"}
    assert ["BAD1"] in calls                                # the bad accession was isolated by recursive halving
    assert calls[0] == ["A1", "A2", "BAD1", "A3"]


def test_fetch_fasta_chunks_requests(monkeypatch):
    calls = []
    monkeypatch.setattr(bf, "_get", lambda b: (calls.append(len(b)) or "".join(f">sp|{a}|N\nM\n" for a in b)))
    seqs = bf.fetch_fasta([f"A{i}" for i in range(7)], chunk=3)
    assert len(seqs) == 7 and calls == [3, 3, 1]


FUNGAMR_HEADER = ["species", "gene or protein name", "mutation", "drug", "confidence score", "ref_seq_uniprot_accession"]


def write_table(fdir, rows):
    fdir.mkdir(parents=True, exist_ok=True)
    with open(fdir / "FungAMR_070425.tsv", "w") as fh:
        fh.write("\t".join(FUNGAMR_HEADER) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def test_main_builds_reference_fasta_and_mutation_table(tmp_path, monkeypatch):
    fdir = tmp_path / "fungamr"
    write_table(fdir, [
        ["Aspergillus fumigatus", "cyp51A", "L98H", "itraconazole", "1", "Q4WNT5"],
        ["Aspergillus fumigatus", "cyp51A", "G54W", "itraconazole", "3", "Q4WNT5"],
        ["Aspergillus fumigatus", "cyp51A", "M220K", "voriconazole", "9", "Q4WNT5"],      # above --min-confidence
        ["Candida albicans", "ERG11", "Y132F", "fluconazole", "2", "P10613, P10613"],    # duplicate accession in one cell
        ["Candida albicans", "FKS1", "S645P", "caspofungin", "-1", "O74166"],           # negative score: not a resistance mutation
        ["Candida glabrata", "PDR1", "", "fluconazole", "nan", "NA"],                     # no accession
    ])
    monkeypatch.setattr(bf, "fetch_fasta", lambda accs, chunk=80: {a: f"SEQ{a}" for a in accs})
    monkeypatch.setattr(sys, "argv", ["build_fungamr_refs.py", "--data-dir", str(tmp_path)])
    bf.main()
    faa = open(fdir / "reference_proteins.faa").read()
    assert ">cyp51A__Q4WNT5__Aspergillus_fumigatus\nSEQQ4WNT5\n" in faa
    assert ">ERG11__P10613__Candida_albicans\n" in faa and ">FKS1__O74166__Candida_albicans\n" in faa
    assert faa.count(">") == 3
    muts = [l.rstrip("\n").split("\t") for l in open(fdir / "fungamr_mutations.tsv")][1:]
    assert [m[3] for m in muts] == ["L98H", "G54W", "Y132F"]          # confidence 1..8 kept; 9 and -1 dropped
    info = json.load(open(fdir / "reference_build.json"))
    assert info == {"accessions_requested": 3, "sequences_written": 3, "mutations": 3, "panel_rows": 2}
