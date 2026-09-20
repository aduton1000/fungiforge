"""Unit tests for bin/extras.py (W2.6) — dbCAN filter, SignalP/EffectorP/PHI-base parsers, mating-type
logic with GenBank synteny, nQuire parsing, and the CLI on canned tool outputs."""
import json, os, subprocess, sys
import pytest
import extras as ex  # noqa: E402  (sys.path set in conftest)


def domtbl_line(target, tlen, query, acc, qlen, evalue, hf, ht):
    f = [target, "-", str(tlen), query, acc, str(qlen), str(evalue), "100.0", "0.1", "1", "1", "1e-10", str(evalue), "50.0", "0.1", str(hf), str(ht), "1", "50", "1", "50", "0.9", "-"]
    return " ".join(f)


def test_parse_domtbl_and_cazyme_filter(tmp_path):
    p = tmp_path / "d.domtbl"
    p.write_text("# comment\n" + "\n".join([
        domtbl_line("P1", 400, "GH18.hmm", "-", 300, 1e-50, 5, 290),       # good
        domtbl_line("P1", 400, "CBM1.hmm", "-", 36, 1e-20, 1, 36),
        domtbl_line("P2", 300, "GT2.hmm", "-", 200, 1e-3, 1, 200),         # weak E-value
        domtbl_line("P3", 300, "AA9.hmm", "-", 230, 1e-30, 1, 40),         # low coverage
        "short line"]) + "\n")
    rows = ex.parse_domtbl(str(p))
    assert len(rows) == 4 and rows[0]["target"] == "P1" and rows[0]["qlen"] == 300
    fams = ex.cazyme_families(rows)
    assert fams == {"P1": ["CBM1", "GH18"]}
    summ = ex.cazyme_summary(fams)
    assert summ["n_cazymes"] == 1 and summ["by_class"] == {"CBM": 1, "GH": 1} and summ["families"] == {"CBM1": 1, "GH18": 1}
    assert ex.parse_domtbl(str(tmp_path / "none")) == []


def test_parse_signalp6_and_effectorp(tmp_path):
    sp = tmp_path / "prediction_results.txt"
    sp.write_text("# SignalP-6.0\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
                  "P1 desc\tSP\t0.01\t0.99\tCS pos: 20-21. Pr: 0.9\nP2\tOTHER\t0.98\t0.02\t\n")
    r = ex.parse_signalp6(str(sp))
    assert r["P1"] == {"prediction": "SP", "signal_peptide": True, "cs": 20} and r["P2"]["signal_peptide"] is False
    ep = tmp_path / "eff.txt"
    ep.write_text("# Identifier\tCytoplasmic effector\tApoplastic effector\tNon-effector\tPrediction\n"
                  "P1 desc\tY (0.8)\t-\t-\tCytoplasmic effector\nP3\t-\t-\tY (0.7)\tNon-effector\nP4      -       Y (0.6)      -       Apoplastic effector\n")
    e = ex.parse_effectorp(str(ep))
    assert e["P1"]["effector"] and not e["P3"]["effector"] and e["P4"] == {"prediction": "Apoplastic effector", "effector": True}
    assert ex.parse_signalp6(None) == {} and ex.parse_effectorp(str(tmp_path / "x")) == {}


def test_parse_phibase_hits_best_hit_and_phenotype(tmp_path):
    p = tmp_path / "phi.tsv"
    p.write_text("\n".join([
        "P1\tQ4WNT5#PHI:1#Cyp51A#746128#Aspergillus_fumigatus#reduced virulence\t60.0\t400\t500\t480\t1e-100\t500\tQ4WNT5#PHI:1#Cyp51A#746128#Aspergillus_fumigatus#reduced virulence",
        "P1\tPHI:2#x#y#unaffected pathogenicity\t80.0\t100\t500\t480\t1e-20\t150\tPHI:2#x#y#unaffected pathogenicity",   # short alignment
        "P2\tPHI:3#gene#sp#loss of pathogenicity\t35.0\t400\t500\t480\t1e-50\t300\tPHI:3#gene#sp#loss of pathogenicity",   # low identity
        "P3\tX1#PHI:4#gene4#1#Sp#increased virulence (hypervirulence)\t55.0\t300\t320\t300\t1e-60\t320\tX1#PHI:4#gene4#1#Sp#increased virulence (hypervirulence)"]) + "\n")
    h = ex.parse_phibase_hits(str(p))
    assert set(h) == {"P1", "P3"} and h["P1"]["phenotype"] == "reduced_virulence" and h["P1"]["gene"] == "Cyp51A" and h["P1"]["organism"] == "Aspergillus_fumigatus"
    assert h["P3"]["gene"] == "gene4" and h["P3"]["organism"] == "Sp"
    assert ex.virulence_summary(h) == {"n_phibase_hits": 2, "by_phenotype": {"increased_virulence_(hypervirulence)": 1, "reduced_virulence": 1}}


def write_gbk(path, genes):
    """genes: list of (locus_tag, start0) on one contig, single CDS each."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    rec = SeqRecord(Seq("A" * 20000), id="c1", name="c1", description="", annotations={"molecule_type": "DNA"})
    for tag, s in genes:
        rec.features.append(SeqFeature(FeatureLocation(s, s + 300, strand=1), type="CDS", qualifiers={"locus_tag": [tag], "protein_id": [f"ncbi:{tag}-T1"]}))
    SeqIO.write(rec, str(path), "genbank")
    return str(path)


def test_mating_type_alpha_hmg_context_and_fallback(tmp_path):
    pytest.importorskip("Bio")
    # window is +-6 genes: FUN_2 sits next to APN2 (FUN_3); FUN_5 is nine genes downstream
    gbk = write_gbk(tmp_path / "g.gbk", [("FUN_1", 100), ("FUN_2", 1000), ("FUN_3", 2000), ("FUN_4", 3000)] + [(f"FILL_{i}", 3500 + 400 * i) for i in range(8)] + [("FUN_5", 15000)])
    rows = [dict(target="FUN_2-T1", acc="PF00505.20", query="HMG_box", evalue=1e-8, tlen=300, qlen=70, hmm_from=1, hmm_to=70, ali_from=1, ali_to=70),
            dict(target="FUN_3-T1", acc="PF03372.24", query="Exo_endo_phos", evalue=1e-30, tlen=300, qlen=200, hmm_from=1, hmm_to=200, ali_from=1, ali_to=200),
            dict(target="FUN_5-T1", acc="PF00505.20", query="HMG_box", evalue=1e-9, tlen=300, qlen=70, hmm_from=1, hmm_to=70, ali_from=1, ali_to=70)]
    m = ex.mating_type(rows, gbk)
    assert m["mating_type"] == "MAT1-2" and m["mat1_2_candidates"] == ["FUN_2-T1"] and m["n_hmg_box_proteins"] == 2   # FUN_5 is outside the window
    alpha = [dict(target="FUN_4-T1", acc="PF04769.14", query="MATalpha_HMGbox", evalue=1e-20, tlen=300, qlen=100, hmm_from=1, hmm_to=100, ali_from=1, ali_to=100)]
    assert ex.mating_type(alpha, gbk)["mating_type"] == "MAT1-1"
    assert ex.mating_type(rows + alpha, gbk)["mating_type"] == "MAT1-1/MAT1-2"
    # no synteny evidence: a single-domain HMG-box protein is accepted as a candidate, others are not
    only_hmg = [rows[2]]
    assert ex.mating_type(only_hmg, gbk)["mating_type"] == "undetermined"
    assert ex.mating_type(only_hmg, None, {"FUN_5-T1": {"PF00505"}})["mating_type"] == "MAT1-2"
    assert ex.mating_type(only_hmg, None, {"FUN_5-T1": {"PF00505", "PF00046"}})["mating_type"] == "undetermined"
    assert ex.mating_type([], gbk)["mating_type"] == "undetermined"


def test_parse_nquire_lrdmodel():
    txt = "file\tfree\tdip\ttri\ttet\td_dip\td_tri\td_tet\nreads.bin\t-1000.5\t-1010.2\t-1200.0\t-1300.1\t9.7\t199.5\t299.6\n"
    r = ex.parse_nquire_lrdmodel(txt)
    assert r["best_model"] == "dip" and r["ploidy_call"] == 2 and r["delta_to_free"] == 9.7
    r2 = ex.parse_nquire_lrdmodel("file\tfree\tdip\ttri\ttet\nx\t-10\t-50\t-12\t-40\n")
    assert r2["ploidy_call"] == 3
    assert ex.parse_nquire_lrdmodel("") is None and ex.parse_nquire_lrdmodel("file\tfree\nx\tbad\n") is None


def test_cli_with_canned_outputs_and_without_tools(tmp_path, helpers):
    prots = helpers["write_fasta"](tmp_path / "p.faa", {"P1-T1": "MKLLA*", "P2-T1": "MSSS", "P3-T1": "MAAA"})
    sp = tmp_path / "sp.txt"; sp.write_text("P1-T1\tSP\t0.01\t0.99\tCS pos: 20-21.\nP2-T1\tOTHER\t0.9\t0.1\t\n")
    ep = tmp_path / "ep.txt"; ep.write_text("P1-T1\tY (0.8)\t-\t-\tCytoplasmic effector\n")
    dt = tmp_path / "d.domtbl"; dt.write_text(domtbl_line("P2-T1", 300, "GH5.hmm", "-", 300, 1e-40, 1, 290) + "\n")
    ph = tmp_path / "ph.tsv"; ph.write_text("P3-T1\tPHI:9#g#s#reduced virulence\t70\t300\t300\t300\t1e-80\t400\tPHI:9#g#s#reduced virulence\n")
    pf = tmp_path / "pf.domtbl"; pf.write_text(domtbl_line("P3-T1", 300, "MATalpha_HMGbox", "-", 100, 1e-20, 1, 100) + "\n")
    out = tmp_path / "e.json"
    cmd = [sys.executable, os.path.join(helpers["BIN"], "extras.py"), "--sample", "S1", "--proteins", prots, "--workdir", str(tmp_path / "w"),
           "--signalp-results", str(sp), "--effectorp-results", str(ep), "--dbcan-domtbl", str(dt), "--phibase-hits", str(ph), "--pfam-domtbl", str(pf), "--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(out))
    assert d["n_proteins"] == 3 and d["n_secreted"] == 1 and d["n_effectors"] == 1 and d["n_cazymes"] == 1 and d["n_phibase_hits"] == 1
    assert d["mating_type"] == "MAT1-1" and d["ploidy"] == "NA" and d["skipped"] == ["nquire"] and "mating=MAT1-1" in r.stdout
    # nothing available: every block null, stage still succeeds and lists what was skipped
    out2 = tmp_path / "e2.json"
    r2 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "extras.py"), "--sample", "S1", "--proteins", prots, "--workdir", str(tmp_path / "w2"), "--out", str(out2)],
                        capture_output=True, text=True, env={**os.environ, "PATH": "/nonexistent"})
    assert r2.returncode == 0, r2.stderr
    d2 = json.load(open(out2))
    assert d2["cazymes"] is None and d2["secretome"] is None and d2["virulence"] is None and d2["mating_type"] == "undetermined"
    assert set(d2["skipped"]) == {"dbcan", "signalp6", "phibase", "mating_hmm", "nquire"}


# ── nQuire ploidy confidence ─────────────────────────────────────────────────
# A real CEA10 run gave 2,491 biallelic sites over 30 Mb and nQuire's best fit was tetraploid,
# for a haploid organism whose k-mer profile said haploid. Density, not fit, decides.

class _FakeRun:
    """Stands in for extras.run(): nQuire create/denoise/view/lrdmodel."""
    def __init__(self, n_sites, tmp):
        self.n_sites, self.tmp = n_sites, tmp

    def __call__(self, cmd, log=None):
        import subprocess
        sub = cmd[1] if len(cmd) > 1 else ""
        if sub == "create":
            open(cmd[cmd.index("-o") + 1] + ".bin", "w").write("x")
        elif sub == "denoise":
            open(cmd[cmd.index("-o") + 1] + ".bin", "w").write("x")
        out = ""
        if sub == "view":
            out = "\n".join("site%d" % i for i in range(self.n_sites))
        elif sub == "lrdmodel":
            out = ("file\tfree\tdip\ttri\ttet\td_dip\td_tri\td_tet\n"
                   "x.bin\t4066.7\t-731.8\t190.1\t3181.0\t4798.5\t3876.6\t885.7\n")
        return subprocess.CompletedProcess(cmd, 0, out, "")


def _ploidy(monkeypatch, tmp_path, n_sites, min_sites=10000):
    monkeypatch.setattr(ex, "run", _FakeRun(n_sites, tmp_path))
    monkeypatch.setattr(ex.shutil, "which", lambda _x: "/usr/bin/nQuire")
    bam = tmp_path / "r.bam"; bam.write_text("bam")
    return ex.ploidy_from_nquire(str(bam), str(tmp_path / "nq"), min_sites=min_sites)


def test_sparse_sites_do_not_yield_a_polyploid_call(monkeypatch, tmp_path):
    r = _ploidy(monkeypatch, tmp_path, 2491)
    assert r["ploidy_call"] == 1 and r["ploidy_confidence"] == "low"
    assert r["model_ploidy_call"] == 4 and r["model_best"] == "tet"   # raw fit preserved
    assert "2491" in r["note"] and "tet" in r["note"]


def test_dense_sites_report_the_model_call(monkeypatch, tmp_path):
    r = _ploidy(monkeypatch, tmp_path, 50000)
    assert r["ploidy_call"] == 4 and r["best_model"] == "tet" and r["ploidy_confidence"] == "ok"


def test_the_threshold_is_configurable(monkeypatch, tmp_path):
    assert _ploidy(monkeypatch, tmp_path, 2491, min_sites=1000)["ploidy_call"] == 4
