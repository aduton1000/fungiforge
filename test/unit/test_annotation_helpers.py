"""Unit tests for bin/annotation_training.py and bin/annotate_stats.py (W2.5)."""
import json, os, subprocess, sys
import annotation_training as at  # noqa: E402  (sys.path set in conftest)
import annotate_stats as ast_  # noqa: E402

MAP = os.path.join(os.path.dirname(__file__), "..", "..", "fungiforge", "resources", "annotation_training.tsv")


def test_bundled_map_is_well_formed():
    rows = at.read_map(MAP)
    assert rows["Aspergillus fumigatus"] == ("aspergillus_fumigatus", "eurotiomycetes") and rows["Candida"] == ("candida_albicans", "saccharomycetes")
    assert all(len(v) == 2 and v[0] and v[1] for v in rows.values()) and len(rows) > 60


def test_choose_species_then_genus_then_default(tmp_path):
    rows = at.read_map(MAP)
    assert at.choose("Aspergillus flavus", rows)["augustus_species"] == "aspergillus_oryzae"
    r = at.choose("Aspergillus lentulus", rows)
    assert (r["augustus_species"], r["busco_db"], r["basis"]) == ("aspergillus_fumigatus", "eurotiomycetes", "genus:Aspergillus")
    d = at.choose("Aspergillus sp.", rows)
    assert d["basis"] == "genus:Aspergillus"
    u = at.choose("unknown", rows)
    assert (u["augustus_species"], u["busco_db"], u["basis"]) == ("anidulans", "dikarya", "default")
    # staged-resource checks: mapped choices are used only when present
    (tmp_path / "aug" / "species" / "anidulans").mkdir(parents=True)
    (tmp_path / "fdb" / "dikarya").mkdir(parents=True)
    c = at.choose("Aspergillus fumigatus", rows, str(tmp_path / "aug"), str(tmp_path / "fdb"))
    assert c["augustus_species"] == "anidulans" and c["busco_db"] == "dikarya" and len(c["notes"]) == 2
    (tmp_path / "aug" / "species" / "aspergillus_fumigatus").mkdir()
    (tmp_path / "fdb" / "eurotiomycetes").mkdir()
    c2 = at.choose("Aspergillus fumigatus", rows, str(tmp_path / "aug"), str(tmp_path / "fdb"))
    assert c2["augustus_species"] == "aspergillus_fumigatus" and c2["busco_db"] == "eurotiomycetes" and c2["notes"] == []


def test_training_cli_prints_shell_assignments(tmp_path, helpers):
    sp = tmp_path / "s.txt"; sp.write_text("Candida auris\thigh\n")
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "annotation_training.py"), "--species", str(sp), "--map", MAP, "--json", str(tmp_path / "t.json")],
                       capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout == "AUGUSTUS_SPECIES=candida_albicans\nBUSCO_DB=saccharomycetes\n"
    assert json.load(open(tmp_path / "t.json"))["basis"] == "genus:Candida"


ANN_HDR = ["GeneID", "TranscriptID", "Feature", "Contig", "Start", "Stop", "Strand", "Name", "Product", "Alias", "EC_number", "BUSCO",
           "PFAM", "InterPro", "EggNog", "COG", "GO Terms", "Secreted", "Membrane", "Protease", "CAZyme", "Notes", "gDNA", "mRNA", "CDS-transcript", "Translation"]


def ann_row(gid, product="hypothetical protein", pfam="", ipr="", egg="", go="", sec="", caz="", prot="", busco="", ec=""):
    d = dict.fromkeys(ANN_HDR, "")
    d.update({"GeneID": gid, "TranscriptID": gid + "-T1", "Feature": "CDS", "Product": product, "PFAM": pfam, "InterPro": ipr, "EggNog": egg,
              "GO Terms": go, "Secreted": sec, "CAZyme": caz, "Protease": prot, "BUSCO": busco, "EC_number": ec})
    return "\t".join(d[k] for k in ANN_HDR)


def test_annotation_stats_and_cli(tmp_path, helpers):
    ann = tmp_path / "S1.annotations.txt"
    ann.write_text("\t".join(ANN_HDR) + "\n" + "\n".join([
        ann_row("FUN_1", "cytochrome P450", pfam="PF00067", ipr="IPR001128", egg="ENOG5:X", go="GO:0005506", sec="SignalP", busco="EOG1"),
        ann_row("FUN_2", pfam="PF00069", go="GO:0004672", caz="GH18"),
        ann_row("FUN_3"),
        ann_row("FUN_4", "Hypothetical protein", prot="S8", ec="3.4.21.-")]) + "\n")
    st = ast_.annotation_stats(str(ann))
    assert st == {"n_annotated": 4, "pct_pfam": 50.0, "pct_interpro": 25.0, "pct_go": 50.0, "pct_eggnog": 25.0, "pct_named_product": 25.0,
                  "n_secreted": 1, "n_cazyme": 1, "n_protease": 1, "n_busco": 1, "n_ec": 1}
    assert ast_.annotation_stats(str(tmp_path / "none.txt")) == {}
    faa = helpers["write_fasta"](tmp_path / "p.faa", {"FUN_1-T1": "MAAA", "FUN_2-T1": "MCCC", "FUN_3-T1": "MDDD", "FUN_4-T1": "MEEE"})
    tr = tmp_path / "training.json"; tr.write_text(json.dumps({"augustus_species": "aspergillus_fumigatus", "busco_db": "eurotiomycetes", "basis": "species:x", "notes": []}))
    out = tmp_path / "a.json"
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "annotate_stats.py"), "--sample", "S1", "--annotations", str(ann), "--proteins", faa,
                        "--training", str(tr), "--eggnog", "yes", "--interproscan", "no", "--genemark", "yes", "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    doc = json.load(open(out))
    assert doc["stage"] == "annotate" and doc["n_proteins"] == 4 and doc["pct_pfam"] == 50.0 and doc["eggnog"] == "yes" and doc["genemark"] == "yes"
    assert doc["training"]["busco_db"] == "eurotiomycetes" and "4 proteins" in r.stdout
