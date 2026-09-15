"""Unit tests for bin/af_resistance.py — panel parsing, ortholog search, substitution / LoF calls,
confidence by polish mode, species relevance, and the cyp51A TR integration (synthetic data)."""
import json, os, subprocess, sys
import pytest

pytest.importorskip("Bio", reason="Biopython is a runtime dependency of the resistance stage (dev extra)")
import af_resistance as afr  # noqa: E402  (sys.path set in conftest)
import cyp51a_TR  # noqa: E402

PANEL = """# synthetic panel
gene\torganism_regex\tdrug_class\tdrugs\tmechanism\thotspot_aa\tknown_mutations\tnote\tsource
cyp51A\tAspergillus fumigatus\tazole\titraconazole,voriconazole\tsubstitution\t54,98,121\tG54W,L98H,Y121F\t-\ttest
cyp51A_promoter\tAspergillus fumigatus\tazole\titraconazole\tpromoter_TR\tNA\tTR34,TR46\t-\ttest
FUR1\tAspergillus spp.\t5-FC\tflucytosine\tloss_of_function\tNA\tNA\t-\ttest
FKS1\tCandida spp.\techinocandin\tcaspofungin\tsubstitution\t639,645\tS639F,S645P\t-\ttest
"""


def make_ref(h):
    """120-aa synthetic cyp51A reference with G at 54, L at 98, Y at 121 (1-based)."""
    p = list(h["random_protein"](130, seed=41))
    p[53], p[97], p[120] = "G", "L", "Y"
    return "".join(p)


def setup(tmp_path, h, ortholog=None, species="Aspergillus fumigatus", extra_prots=None, ref=True, fur1=None):
    panel = tmp_path / "panel.tsv"; panel.write_text(PANEL)
    cyp = make_ref(h)
    fur = h["random_protein"](200, seed=42)
    refs = {"cyp51A__Q4WNT5__Aspergillus_fumigatus": cyp, "FUR1__Q9UVL5__Aspergillus_fumigatus": fur,
            "Cdr1_Erg11_Fcy1_Fks1__P0ABC__Candida_albicans": h["random_protein"](150, seed=43)}
    ref_faa = h["write_fasta"](tmp_path / "ref.faa", refs) if ref else None
    prots = dict(extra_prots or {})
    if ortholog is not None:
        prots["FUN_000123-T1 cyp51A-like"] = ortholog
    if fur1 is not None:
        prots["FUN_000200-T1"] = fur1
    prots["FUN_000999-T1"] = h["random_protein"](300, seed=44)          # unrelated protein
    proteins = h["write_fasta"](tmp_path / "prot.faa", prots)
    sp = tmp_path / "species.txt"; sp.write_text(f"{species}\thigh\n")
    return {"panel": str(panel), "ref": ref_faa, "proteins": proteins, "species": str(sp), "cyp": cyp, "fur": fur}


def run(h, tmp_path, s, polish="hybrid", assembly=None, gbk=None, data_dir=""):
    out = tmp_path / "res.json"
    cmd = [sys.executable, os.path.join(h["BIN"], "af_resistance.py"), "--sample", "S1", "--proteins", s["proteins"],
           "--species", s["species"], "--panel", s["panel"], "--polish-mode", polish, "--out", str(out), "--data-dir", data_dir]
    if s["ref"]: cmd += ["--ref-faa", s["ref"]]
    if assembly: cmd += ["--assembly", assembly]
    if gbk: cmd += ["--gbk", gbk]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.load(open(out)), r


def by_gene(doc, gene):
    return [c for c in doc["calls"] if c["gene"] == gene]


def sub(tmp_path, name):
    d = tmp_path / name; d.mkdir(); return d


# ---------------------------------------------------------------- helpers
def test_read_panel_skips_comments_and_maps_columns(tmp_path):
    p = tmp_path / "p.tsv"; p.write_text(PANEL)
    rows = afr.read_panel(p)
    assert [r["gene"] for r in rows] == ["cyp51A", "cyp51A_promoter", "FUR1", "FKS1"]
    assert rows[0]["known_mutations"] == "G54W,L98H,Y121F" and rows[2]["mechanism"] == "loss_of_function"


def test_reference_index_indexes_every_gene_keyword_in_a_header(tmp_path, helpers):
    faa = helpers["write_fasta"](tmp_path / "r.faa", {"Cdr1_Erg11_Fcy1_Fks1__ACC__Candida_albicans": "MKV", "cyp51A__X__Af": "MAA", "CYP51__Y__Ca": "MCC"})
    idx = afr.load_reference_index("", faa, None)
    assert set(idx) >= {"cdr1", "erg11", "fcy1", "fks1", "cyp51a", "cyp51"}
    assert idx["cyp51a"][0][2] == "MAA" and idx["cyp51"][0][2] == "MCC"


def test_best_ortholog_and_residue_mapping(helpers):
    ref = make_ref(helpers)
    q = ref[:97] + "H" + ref[98:]
    prots = {"orth": q, "junk": helpers["random_protein"](130, seed=99)}
    name, seq, pid, aln = afr.best_ortholog(ref, prots)
    assert name == "orth" and pid > 0.99
    assert afr.ref_to_query_residue(aln, 98) == "H" and afr.ref_to_query_residue(aln, 54) == "G"


def test_best_ortholog_returns_none_below_identity_floor(helpers):
    ref = make_ref(helpers)
    assert afr.best_ortholog(ref, {"x": helpers["random_protein"](130, seed=77)}) is None


# ---------------------------------------------------------------- end to end
def test_known_mutation_is_called_with_high_confidence_for_hybrid(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref[:97] + "H" + ref[98:])
    doc, r = run(helpers, tmp_path, s, polish="hybrid")
    c = by_gene(doc, "cyp51A")
    assert len(c) == 1 and c[0]["status"] == "variant" and c[0]["change"] == "L98H" and c[0]["known"] is True
    assert c[0]["class"] == "known_resistance_mutation" and c[0]["confidence"] == "high"
    assert doc["summary"]["resistant_drug_classes"] == ["azole"] and doc["summary"]["n_known_mutations"] == 1
    assert doc["summary"]["reference_available"] is True and "known=1" in r.stdout


def test_wild_type_ortholog_and_novel_hotspot_variant(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref)
    doc, _ = run(helpers, tmp_path, s)
    assert by_gene(doc, "cyp51A")[0]["status"] == "wild_type" and doc["summary"]["resistant_drug_classes"] == []
    b = sub(tmp_path, "b")
    doc2, _ = run(helpers, b, setup(b, helpers, ortholog=ref[:53] + "A" + ref[54:]))
    c = by_gene(doc2, "cyp51A")[0]
    assert c["status"] == "variant" and c["change"] == "G54A" and c["known"] is False and c["class"] == "novel_hotspot_variant"
    assert doc2["summary"]["resistant_drug_classes"] == []                 # novel variants do not claim resistance


def test_ont_only_calls_are_provisional(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref[:97] + "H" + ref[98:])
    doc, _ = run(helpers, tmp_path, s, polish="ont_only")
    assert by_gene(doc, "cyp51A")[0]["confidence"] == "provisional_ont_only" and doc["polish_mode"] == "ont_only"


def test_not_detected_when_no_ortholog(tmp_path, helpers):
    s = setup(tmp_path, helpers, ortholog=None)
    doc, _ = run(helpers, tmp_path, s)
    assert by_gene(doc, "cyp51A")[0]["status"] == "not_detected"


def test_no_reference_is_reported_not_guessed(tmp_path, helpers):
    ref = make_ref(helpers)
    s = setup(tmp_path, helpers, ortholog=ref, ref=False)
    doc, _ = run(helpers, tmp_path, s)
    c = by_gene(doc, "cyp51A")[0]
    assert c["status"] == "no_reference" and "install FungAMR" in c["note"]
    assert doc["summary"]["reference_available"] is False


def test_loss_of_function_by_truncation(tmp_path, helpers):
    s = setup(tmp_path, helpers, ortholog=make_ref(helpers))
    fur = s["fur"]
    t = sub(tmp_path, "t")
    doc, _ = run(helpers, t, setup(t, helpers, ortholog=make_ref(helpers), fur1=fur[:90]))
    c = by_gene(doc, "FUR1")[0]
    assert c["status"] == "loss_of_function" and "truncated" in c["note"]
    assert "5-FC" in doc["summary"]["resistant_drug_classes"]
    f = sub(tmp_path, "f")
    doc2, _ = run(helpers, f, setup(f, helpers, ortholog=make_ref(helpers), fur1=fur))
    assert by_gene(doc2, "FUR1")[0]["status"] == "intact"


def test_species_relevance_filters_panel_rows(tmp_path, helpers):
    s = setup(tmp_path, helpers, ortholog=make_ref(helpers))
    doc, _ = run(helpers, tmp_path, s)
    assert doc["genes_searched"] == ["FUR1", "cyp51A", "cyp51A_promoter"]      # Candida FKS1 row not searched
    c = sub(tmp_path, "c")
    doc2, _ = run(helpers, c, setup(c, helpers, species="Candida auris"))
    assert doc2["genes_searched"] == ["FKS1"] and doc2["cyp51A_TR"] is None
    u = sub(tmp_path, "u")
    doc3, _ = run(helpers, u, setup(u, helpers, species="unknown"))
    assert doc3["genes_searched"] == [] and doc3["calls"] == []


def test_tr34_detected_via_assembly_and_gbk(tmp_path, helpers):
    unit = helpers["TR34_UNIT"]
    promoter = helpers["random_dna"](532, seed=51) + unit * 2
    gene = helpers["random_dna"](900, seed=52)
    seq = helpers["random_dna"](100, seed=53) + promoter + gene + helpers["random_dna"](100, seed=54)
    gs = 100 + len(promoter)
    gbk = helpers["write_gbk"](tmp_path / "g.gbk", seq, [(gs, gs + len(gene), 1, "gene", {"gene": ["cyp51A"]})])
    asm = helpers["write_fasta"](tmp_path / "a.fa", {"contig_1": seq})
    s = setup(tmp_path, helpers, ortholog=make_ref(helpers))
    doc, _ = run(helpers, tmp_path, s, polish="hybrid", assembly=asm, gbk=gbk)
    tr = by_gene(doc, "cyp51A_promoter")
    assert len(tr) == 1 and tr[0]["status"] == "resistance" and tr[0]["change"] == "TR34" and tr[0]["confidence"] == "high"
    assert doc["cyp51A_TR"]["tr_type"] == "TR34" and doc["summary"]["resistant_drug_classes"] == ["azole"]
    doc_ont, _ = run(helpers, tmp_path, s, polish="ont_only", assembly=asm, gbk=gbk)
    assert by_gene(doc_ont, "cyp51A_promoter")[0]["confidence"] == "medium"


def test_missing_panel_fails_loudly(tmp_path, helpers):
    s = setup(tmp_path, helpers, ortholog=make_ref(helpers))
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "af_resistance.py"), "--sample", "S1", "--proteins", s["proteins"],
                        "--species", s["species"], "--panel", str(tmp_path / "missing.tsv"), "--out", str(tmp_path / "o.json")],
                       capture_output=True, text=True)
    # the packaged panel is a legitimate fallback when the fungiforge package is importable; otherwise it must exit non-zero
    assert r.returncode != 0 or "using packaged copy" in r.stderr
