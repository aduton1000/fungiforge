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


# funannotate 1.8.17 writes "mRNA" in the Feature column, not "CDS". The fixture said CDS, which
# is why every test passed while a real 9,844-gene table produced n_annotated = 0.
def ann_row(gid, product="hypothetical protein", pfam="", ipr="", egg="", go="", sec="", caz="", prot="", busco="", ec="",
            feature="mRNA"):
    d = dict.fromkeys(ANN_HDR, "")
    d.update({"GeneID": gid, "TranscriptID": gid + "-T1", "Feature": feature, "Product": product, "PFAM": pfam, "InterPro": ipr, "EggNog": egg,
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
    expect = {"n_annotated": 4, "pct_pfam": 50.0, "pct_interpro": 25.0, "pct_go": 50.0, "pct_eggnog": 25.0,
              "pct_named_product": 25.0, "n_secreted": 1, "n_cazyme": 1, "n_protease": 1, "n_busco": 1, "n_ec": 1}
    assert {k: st[k] for k in expect} == expect
    assert st["annotations_header"] == ANN_HDR and "annotations_note" not in st
    missing = ast_.annotation_stats(str(tmp_path / "none.txt"))
    assert "pct_pfam" not in missing and "not found" in missing["annotations_note"]
    faa = helpers["write_fasta"](tmp_path / "p.faa", {"FUN_1-T1": "MAAA", "FUN_2-T1": "MCCC", "FUN_3-T1": "MDDD", "FUN_4-T1": "MEEE"})
    tr = tmp_path / "training.json"; tr.write_text(json.dumps({"augustus_species": "aspergillus_fumigatus", "busco_db": "eurotiomycetes", "basis": "species:x", "notes": []}))
    out = tmp_path / "a.json"
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "annotate_stats.py"), "--sample", "S1", "--annotations", str(ann), "--proteins", faa,
                        "--training", str(tr), "--eggnog", "yes", "--interproscan", "no", "--genemark", "yes", "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    doc = json.load(open(out))
    assert doc["stage"] == "annotate" and doc["n_proteins"] == 4 and doc["pct_pfam"] == 50.0 and doc["eggnog"] == "yes" and doc["genemark"] == "yes"
    assert doc["training"]["busco_db"] == "eurotiomycetes" and "4 proteins" in r.stdout


# ── why the coverage columns are empty ───────────────────────────────────────
# All four pct_* columns came back NA on a real run while eggNOG and InterProScan had both
# done real work. annotation_stats returned a bare {}, so the master row could not say whether
# the table was missing, empty, or filtered away. Every empty result now carries the reason.

def test_missing_annotations_table_explains_itself():
    st = ast_.annotation_stats("/no/such/file.annotations.txt")
    assert "pct_pfam" not in st
    assert "not found" in st["annotations_note"] and "file.annotations.txt" in st["annotations_note"]


def test_no_annotations_argument_explains_itself():
    assert "annotations_note" in ast_.annotation_stats(None)


def test_header_only_table_reports_zero_and_the_header(tmp_path):
    p = tmp_path / "S.annotations.txt"
    p.write_text("\t".join(ANN_HDR) + "\n")
    st = ast_.annotation_stats(str(p))
    assert st["n_annotated"] == 0 and st["annotations_header"] == ANN_HDR
    assert "no usable rows" in st["annotations_note"]


def test_table_without_a_feature_column_is_still_counted(tmp_path):
    """A different funannotate layout must not filter every row away and report nothing."""
    hdr = [c for c in ANN_HDR if c != "Feature"]
    p = tmp_path / "S.annotations.txt"
    row = dict(zip(ANN_HDR, ann_row("FUN_1", "cytochrome P450", pfam="PF00067").split("\t")))
    p.write_text("\t".join(hdr) + "\n" + "\t".join(row[c] for c in hdr) + "\n")
    st = ast_.annotation_stats(str(p))
    assert st["n_annotated"] == 1 and st["pct_pfam"] == 100.0


# ── which Feature values count as protein-coding ─────────────────────────────

def _table(tmp_path, rows):
    p = tmp_path / "S.annotations.txt"
    p.write_text("\t".join(ANN_HDR) + "\n" + "\n".join(rows) + "\n")
    return str(p)


def test_mrna_and_cds_rows_both_count(tmp_path):
    """The label differs between funannotate builds; both mean a protein-coding transcript."""
    for feature in ("mRNA", "CDS", ""):
        st = ast_.annotation_stats(_table(tmp_path, [
            ann_row("G1", "cytochrome P450", pfam="PF00067", feature=feature),
            ann_row("G2", feature=feature)]))
        assert st["n_annotated"] == 2, feature
        assert st["pct_pfam"] == 50.0


def test_noncoding_rows_are_excluded(tmp_path):
    st = ast_.annotation_stats(_table(tmp_path, [
        ann_row("G1", "cytochrome P450", pfam="PF00067", feature="mRNA"),
        ann_row("T1", feature="tRNA"),
        ann_row("R1", feature="rRNA")]))
    assert st["n_annotated"] == 1 and st["pct_pfam"] == 100.0


def test_an_unfamiliar_feature_label_still_counts(tmp_path):
    """Better to count an unknown label than to silently report nothing, as CDS-only filtering did."""
    st = ast_.annotation_stats(_table(tmp_path, [ann_row("G1", pfam="PF00067", feature="transcript")]))
    assert st["n_annotated"] == 1
