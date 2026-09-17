"""Unit tests for bin/fungamr_panel.py — panel rows and the evidence map derived from FungAMR."""
import os
import fungamr_panel as fp  # noqa: E402  (sys.path set in conftest)

HEADER = ["species", "gene or protein name", "mutation", "drug", "confidence score", "pubmedid"]


def rows(*data):
    return [dict(zip(HEADER, d)) for d in data]


def test_derive_panel_tiers_classes_synonyms_and_filters():
    panel, ev = fp.derive_panel(rows(
        ["Aspergillus fumigatus", "Cyp51A", "L98H,TR34", "Itraconazole", "1", "111"],
        ["Aspergillus fumigatus", "Cyp51A", "L98H", "Voriconazole", "8", "222"],        # weaker tier for the same change: keep 1
        ["Aspergillus fumigatus", "Cyp51A", "G54W", "Itraconazole", "8", "333"],
        ["Aspergillus fumigatus", "Cyp51A", "S297T", "Itraconazole", "-3", "444"],      # negative: shown not to confer
        ["Aspergillus fumigatus", "Cyp51A", "Deletion", "Itraconazole", "5", "555"],    # not a substitution
        ["Aspergillus fumigatus", "Cyp51A", "M220I", "Itraconazole", "nan", "666"],     # no score
        ["Candidozyma auris", "Erg11", "Y132F", "Fluconazole", "2", "777"],
        ["Candidozyma auris", "Fks1", "S639F", "Caspofungin", "8", "888"],
        ["Candidozyma auris", "Fks1", "S639F", "Micafungin", "4", "999"],
    ), fp.DEFAULT_CLASS)
    assert [(r["gene"], r["organism"]) for r in panel] == [("Cyp51A", "Aspergillus fumigatus"), ("Erg11", "Candidozyma auris"), ("Fks1", "Candidozyma auris")]
    cyp = panel[0]
    assert cyp["known_mutations"] == "G54W,L98H" and cyp["hotspot_aa"] == "54,98" and cyp["tiers"] == "G54W:8,L98H:1" and cyp["mechanism"] == "substitution"
    assert cyp["drug_class"] == "azole" and "111" in cyp["source"]
    assert ev[("cyp51a", "aspergillus fumigatus", "L98H")] == {"tier": 1, "classes": ["azole"]}
    assert ev[("erg11", "candida auris", "Y132F")]["tier"] == 2 and ev[("fks1", "candidozyma auris", "S639F")] == {"tier": 4, "classes": ["echinocandin"]}
    assert panel[1]["organism_regex"] == r"Candidozyma\ auris|Candida\ auris"
    panel7, _ = fp.derive_panel(rows(["Aspergillus fumigatus", "Cyp51A", "G54W", "Itraconazole", "8", "1"]), fp.DEFAULT_CLASS, min_tier=7)
    assert panel7 == []


def test_drug_classes_from_csv_and_defaults(tmp_path):
    p = tmp_path / "drugs_class.csv"
    p.write_text("drug;usage;class\nTebuconazole;Agricultural;Agricultural azoles\nRezafungin;Clinical;Echinocandins\nOlorofim;Clinical;Orotomides\n")
    c = fp.drug_classes(str(p))
    assert c["tebuconazole"] == "azole_agricultural" and c["rezafungin"] == "echinocandin" and c["olorofim"] == "orotomide" and c["fluconazole"] == "azole"
    assert fp.drug_classes(str(tmp_path / "none.csv"))["caspofungin"] == "echinocandin"


def test_write_then_load_roundtrip_and_on_the_fly_derivation(tmp_path):
    fdir = tmp_path / "fungamr"; fdir.mkdir()
    (fdir / "FungAMR_070425.tsv").write_text("\t".join(HEADER) + "\n" + "\n".join("\t".join(r) for r in [
        ["Candida albicans", "Erg11", "Y132H", "Fluconazole", "2", "1"], ["Candida albicans", "Erg11", "K143R", "Fluconazole", "8", "2"]]) + "\n")
    panel, ev, src = fp.load_or_derive(str(tmp_path))
    assert src == "FungAMR_070425.tsv" and len(panel) == 1 and ev[("erg11", "candida albicans", "K143R")]["tier"] == 8
    fp.write_panel(panel, str(fdir / "fungamr_panel.tsv"))
    panel2, ev2, src2 = fp.load_or_derive(str(tmp_path))
    assert src2 == "fungamr_panel.tsv" and panel2[0]["known_mutations"] == "Y132H,K143R" and ev2 == ev
    assert fp.load_or_derive("") == ([], {}, None) and fp.load_or_derive(str(tmp_path / "nowhere")) == ([], {}, None)


def test_cli_writes_panel(tmp_path, helpers):
    import subprocess, sys
    fdir = tmp_path / "fungamr"; fdir.mkdir()
    (fdir / "FungAMR_1.tsv").write_text("\t".join(HEADER) + "\nCandida albicans\tErg11\tY132H\tFluconazole\t2\t1\n")
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "fungamr_panel.py"), "--data-dir", str(tmp_path)], capture_output=True, text=True)
    assert r.returncode == 0 and "1 species/gene rows" in r.stdout and (fdir / "fungamr_panel.tsv").exists()
