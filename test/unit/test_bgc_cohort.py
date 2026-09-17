"""Unit tests for the W3.2 pieces: KnownClusterBlast parsing + mycotoxin flags in bin/bgc_summary.py,
and bin/bgc_families.py (region GenBank parsing, DIAMOND similarity, family clustering, CLI)."""
import json, os, subprocess, sys
import bgc_summary as bg  # noqa: E402  (sys.path set in conftest)
import bgc_families as bf  # noqa: E402

COMP = os.path.join(os.path.dirname(__file__), "..", "..", "fungiforge", "resources", "mycotoxin_compounds.tsv")


def as_json(tmp_path, kcb_key="knowncluster"):
    rec = {"id": "contig_1", "features": [
        {"type": "region", "location": "[1000:60000]", "qualifiers": {"region_number": ["1"], "product": ["NRPS"]}},
        {"type": "region", "location": "[90000:130000]", "qualifiers": {"region_number": ["2"], "product": ["T1PKS"]}}],
        "modules": {"antismash.modules.clusterblast": {kcb_key: {"results": [
            {"region_number": 1, "total_hits": 3, "ranking": [[{"accession": "BGC0000361", "cluster_label": "gliotoxin", "proteins": ["a", "b", "c", "d"], "description": "gliotoxin", "cluster_type": "NRP", "tags": []},
                                                              {"core_gene_hits": 2, "blast_score": 1500, "synteny_score": 4, "core_bonus": 3, "pairings": [["q1", 0, {}], ["q2", 1, {}]]}]]},
            {"region_number": 2, "total_hits": 0, "ranking": []}]}}}}
    d = tmp_path / "as"; d.mkdir(parents=True, exist_ok=True)
    (d / "S1.json").write_text(json.dumps({"records": [rec]}))
    return str(d)


def test_knowncluster_hits_and_flags(tmp_path):
    hits = bg.knowncluster_hits(as_json(tmp_path))
    assert set(hits) == {"contig_1|1", "contig_1|2"} and hits["contig_1|2"]["best"] is None
    b = hits["contig_1|1"]["best"]
    assert b["accession"] == "BGC0000361" and b["description"] == "gliotoxin" and b["n_query_proteins_hit"] == 2 and b["n_reference_proteins"] == 4
    flags = bg.flag_compounds(hits, bg.load_compounds(COMP))
    assert len(flags) == 1 and flags[0]["compound"] == "gliotoxin" and flags[0]["class"] == "mycotoxin" and flags[0]["region"] == "contig_1|1"
    assert bg.knowncluster_hits(as_json(tmp_path / "alt", kcb_key="knownclusterblast"))["contig_1|1"]["best"]["accession"] == "BGC0000361"
    comps = bg.load_compounds(COMP)
    assert any(n == "aflatoxin" for n, _, _ in comps) and all(c in ("mycotoxin", "antifungal_antibiotic", "other_bioactive") for _, c, _ in comps)


def test_bgc_summary_cli_reports_flags(tmp_path, helpers):
    d = as_json(tmp_path)
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "bgc_summary.py"), "--sample", "S1", "--as-dir", d, "--status", "ok", "--out", str(tmp_path / "b.json")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    doc = json.load(open(tmp_path / "b.json"))
    assert doc["n_clusters"] == 2 and doc["knownclusterblast"]["ran"] and doc["knownclusterblast"]["n_regions_with_mibig_hit"] == 1
    assert doc["mycotoxin_compounds"] == ["gliotoxin"] and doc["bioactive_compounds"] == [] and "mycotoxins gliotoxin" in r.stdout
    (tmp_path / "as2").mkdir(); (tmp_path / "as2" / "S1.json").write_text(json.dumps({"records": [{"id": "c", "features": [{"type": "region", "qualifiers": {"region_number": ["1"], "product": ["NRPS"]}}], "modules": {}}]}))
    r2 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "bgc_summary.py"), "--sample", "S1", "--as-dir", str(tmp_path / "as2"), "--status", "ok", "--out", str(tmp_path / "b2.json")], capture_output=True, text=True)
    d2 = json.load(open(tmp_path / "b2.json"))
    assert r2.returncode == 0 and d2["knownclusterblast"]["ran"] is False and "cb-knownclusters" in d2["knownclusterblast"]["note"] and d2["mycotoxin_compounds"] == []


def region_gbk(records):
    """records: [(name, region_number, product, {cds_name: protein})] -> GenBank text."""
    out = []
    for name, rn, prod, cds in records:
        out.append(f"LOCUS       {name}  60000 bp    DNA     linear   PLN 01-JAN-2026\nFEATURES             Location/Qualifiers\n")
        out.append(f"     region          1..60000\n                     /region_number=\"{rn}\"\n                     /product=\"{prod}\"\n")
        for cn, seq in cds.items():
            out.append(f"     CDS             1..300\n                     /locus_tag=\"{cn}\"\n                     /translation=\"{seq[:40]}\n                     {seq[40:]}\"\n")
        out.append("ORIGIN\n        1 acgt\n//\n")
    return "".join(out)


def test_parse_regions_similarity_and_families(tmp_path, helpers):
    p = tmp_path / "S1.regions.gbk"
    p.write_text(region_gbk([("c1", 1, "NRPS", {"g1": "M" + "A" * 60, "g2": "M" + "C" * 60, "g3": "M" + "D" * 60}), ("c2", 2, "T1PKS", {"h1": "M" + "E" * 60})]))
    q = tmp_path / "S2.regions.gbk"
    q.write_text(region_gbk([("c9", 1, "NRPS", {"k1": "M" + "A" * 60, "k2": "M" + "C" * 60, "k3": "M" + "F" * 60, "k4": "M" + "G" * 60})]))
    regs = bf.parse_regions(str(p), "S1") + bf.parse_regions(str(q), "S2")
    assert [r["id"] for r in regs] == ["S1|c1|r1", "S1|c2|r2", "S2|c9|r1"] and regs[0]["product"] == ["NRPS"] and len(regs[0]["proteins"]) == 3 and regs[0]["proteins"]["g1"] == "M" + "A" * 60
    dm = "\n".join(["S1|c1|r1||g1\tS2|c9|r1||k1\t100\t61\t61", "S2|c9|r1||k1\tS1|c1|r1||g1\t100\t61\t61",
                    "S1|c1|r1||g2\tS2|c9|r1||k2\t95\t61\t61", "S2|c9|r1||k2\tS1|c1|r1||g2\t95\t61\t61",
                    "S1|c1|r1||g1\tS1|c1|r1||g1\t100\t61\t61", "S1|c2|r2||h1\tS2|c9|r1||k4\t40\t61\t61"]) + "\n"   # self hit ignored; 40 % identity ignored
    hits = bf.parse_diamond(dm)
    sim = bf.region_similarity(regs, hits)
    assert sim == {(0, 2): 2 / 3}                        # 2 of the smaller region's 3 proteins shared
    assert bf.families(regs, sim, 0.5) == [[0, 2], [1]] and bf.families(regs, sim, 0.7) == [[0], [1], [2]]
    dmf = tmp_path / "dm.tsv"; dmf.write_text(dm)
    bj = tmp_path / "S1.bgc.json"; bj.write_text(json.dumps({"knownclusterblast": {"best_hits": {"c1|1": {"accession": "BGC0000361", "description": "gliotoxin", "core_gene_hits": 2, "blast_score": 900}}}}))
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "bgc_families.py"), "--outdir", str(tmp_path / "cohort"), "--regions", str(p), str(q), "--bgc", str(bj), "--diamond-hits", str(dmf)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(tmp_path / "cohort" / "bgc_families.json"))
    assert d["n_regions"] == 3 and d["n_families"] == 2 and d["n_families_shared"] == 1 and d["known_families"] == ["GCF0001"]
    f = d["families"][0]
    assert f["family"] == "GCF0001" and f["isolates"] == ["S1", "S2"] and f["mibig_best"]["description"] == "gliotoxin" and f["products"] == ["NRPS"]
    tsv = (tmp_path / "cohort" / "bgc_families.tsv").read_text().splitlines()
    assert tsv[0] == "family\tregion\tsample\tproducts\tn_proteins\tmibig_best" and tsv[1].startswith("GCF0001\tS1|c1|r1\tS1\tNRPS\t3\tgliotoxin")
    r2 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "bgc_families.py"), "--outdir", str(tmp_path / "c2"), "--regions", str(tmp_path / "missing.gbk")], capture_output=True, text=True)
    assert r2.returncode == 0 and json.load(open(tmp_path / "c2" / "bgc_families.json"))["families"] == []
