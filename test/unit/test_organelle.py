"""Unit tests for bin/mito_extract.py and bin/mito_annotate.py (W2.7)."""
import json, os, subprocess, sys
import mito_extract as mx  # noqa: E402  (sys.path set in conftest)
import mito_annotate as ma  # noqa: E402


def hit(q, contig, sstart, send, ev=1e-30, qs=1, qe=100, bits=200.0, qlen=None, pid=60.0):
    f = [q, contig, str(pid), str(abs(send - sstart) // 3 + 1), str(qs), str(qe), str(sstart), str(send), str(ev), str(bits)]
    if qlen is not None:
        f.append(str(qlen))
    return "\t".join(f)


def test_extract_classifies_mito_contigs_by_genes_length_and_gc(tmp_path, helpers):
    at_rich = ("AT" * 100 + "GC" * 10) * 60                     # 13.2 kb, GC 9 %
    seqs = {"mt": at_rich, "chr1": helpers["random_dna"](60000, seed=3), "numt": helpers["random_dna"](30000, seed=4), "small": at_rich[:5000]}
    fa = helpers["write_fasta"](tmp_path / "a.fa", seqs)
    t = tmp_path / "t.tsv"
    t.write_text("\n".join([hit("cox1|X|Af", "mt", 100, 1600), hit("cob|X|Af", "mt", 3000, 4100), hit("nad5|X|Af", "chr1", 10, 1500), hit("cox1|X|Af", "chr1", 5000, 6000),
                            hit("cox2|X|Af", "numt", 100, 800), hit("cob|X|Af", "numt", 2000, 3000), hit("atp6|X|Af", "small", 10, 700, ev=1e-3)]) + "\n")
    genes = mx.genes_per_contig(str(t))
    assert genes == {"mt": {"cob", "cox1"}, "chr1": {"cox1", "nad5"}, "numt": {"cob", "cox2"}}          # weak hit dropped
    mito, why = mx.classify(mx.read_fasta(fa), genes, max_len=50000, min_genes=2, max_gc=0.40)
    assert mito == ["mt"] and why["chr1"]["reason"] == "too long" and why["numt"]["reason"].startswith("GC too high")
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "mito_extract.py"), "--sample", "S1", "--assembly", fa, "--tblastn", str(t), "--max-len", "50000",
                        "--out-mito", str(tmp_path / "m.fa"), "--out-nuclear", str(tmp_path / "n.fa"), "--json", str(tmp_path / "x.json")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(tmp_path / "x.json"))
    assert d["mito_contigs"] == ["mt"] and d["mito_bp"] == len(at_rich) and d["nuclear_contigs"] == 3
    assert mx.read_fasta(str(tmp_path / "m.fa")) == {"mt": at_rich} and set(mx.read_fasta(str(tmp_path / "n.fa"))) == {"chr1", "numt", "small"}


def test_annotate_genes_chains_exons_counts_introns_and_copies(tmp_path):
    rows = ma.parse_hits(None)
    assert rows == []
    t = tmp_path / "g.tsv"
    t.write_text("\n".join([
        hit("cox1|A|Af293", "mt", 1000, 1600, qs=1, qe=200, bits=300, qlen=530),        # exon 1
        hit("cox1|A|Af293", "mt", 2800, 3800, qs=201, qe=530, bits=500, qlen=530),      # exon 2 after a 1.2 kb intron
        hit("cox1|B|Scer", "mt", 1000, 3800, qs=1, qe=500, bits=400, qlen=520),          # weaker species
        hit("cox1|A|Af293", "mt2", 100, 700, qs=1, qe=200, bits=150, qlen=530),          # second copy elsewhere
        hit("cob|A|Af293", "mt", 12000, 13100, qs=1, qe=385, bits=600, qlen=385),
        hit("atp8|A|Af293", "mt", 5000, 5150, qs=1, qe=48, bits=20, qlen=48, ev=1e-2)]) + "\n")   # weak
    g = ma.annotate_genes(ma.parse_hits(str(t)))
    assert g["cox1"]["present"] and g["cox1"]["reference"] == "Af293" and g["cox1"]["n_exons"] == 2 and g["cox1"]["n_introns"] == 1
    assert g["cox1"]["intron_lengths"] == [1199] and g["cox1"]["query_coverage"] == 1.0 and (g["cox1"]["start"], g["cox1"]["end"]) == (1000, 3800)
    assert g["cox1"]["copies"] == 1 and g["cox1"]["contigs_with_gene"] == 2          # a second contig carries it: fragment, not a duplication
    assert g["cob"]["n_introns"] == 0 and g["cob"]["copies"] == 1 and g["atp8"] == {"present": False} and g["nad1"] == {"present": False}
    r = tmp_path / "r.tsv"
    r.write_text(hit("rnl|X|Af", "mt", 20000, 24500, pid=88.0) + "\n" + hit("rns|X|Af", "mt", 6000, 6100, pid=90.0) + "\n")
    rr = ma.annotate_rrna(ma.parse_hits(str(r), with_qlen=False))
    assert rr["rnl"]["present"] and rr["rnl"]["strand"] == "+" and rr["rns"] == {"present": False}     # rns alignment too short


def test_circular_and_heteroplasmy(tmp_path, helpers):
    core = helpers["random_dna"](3000, seed=8)
    assert ma.circular(core + core[:120]) == 120 and ma.circular(core) == 0 and ma.circular(core + core[:30]) == 0
    p = tmp_path / "m.pileup"
    p.write_text("mt\t1\tA\t30\t" + "." * 30 + "\t" + "I" * 30 + "\n"
                 "mt\t2\tC\t40\t" + "." * 30 + "T" * 10 + "\t" + "I" * 40 + "\n"         # 25 % minor allele
                 "mt\t3\tG\t10\t" + "," * 5 + "a" * 5 + "\tIIIIIIIIII\n"                  # too shallow
                 "mt\t4\tT\t25\t" + "." * 22 + "^]." + "+2AA" + ".$" + "\t" + "I" * 25 + "\n")   # read starts / indel / read end handled
    h = ma.heteroplasmy(str(p))
    assert h["n_positions_covered"] == 3 and h["n_heteroplasmic_sites"] == 1 and h["sites"][0]["pos"] == 2 and h["sites"][0]["alleles"] == {"ref": 0.75, "T": 0.25}
    assert ma.heteroplasmy(str(p), min_frac=0.30)["n_heteroplasmic_sites"] == 0                      # ONT floor
    q = tmp_path / "i.pileup"
    q.write_text("mt\t9\tA\t30\t" + "." * 22 + "+3ACG" * 8 + "\t" + "I" * 30 + "\n")               # indels only: not heteroplasmy
    assert ma.heteroplasmy(str(q))["n_heteroplasmic_sites"] == 0
    assert ma.heteroplasmy(str(tmp_path / "none")) is None


def test_annotate_cli_with_and_without_mito(tmp_path, helpers):
    core = helpers["random_dna"](8000, seed=9)
    fa = helpers["write_fasta"](tmp_path / "m.fa", {"mt": core + core[:80]})
    t = tmp_path / "g.tsv"; t.write_text(hit("cox1|A|Af293", "mt", 100, 1600, qs=1, qe=500, bits=500, qlen=530) + "\n")
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "mito_annotate.py"), "--sample", "S1", "--mito", fa, "--tblastn", str(t), "--mito-depth", "800", "--nuclear-depth", "40",
                        "--out-json", str(tmp_path / "o.json"), "--out-gff", str(tmp_path / "o.gff")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(tmp_path / "o.json"))
    assert d["stage"] == "organelle" and d["mito_present"] and d["mito_size"] == 8080 and d["n_core_genes"] == 1 and d["circular"] is True and d["copy_ratio"] == 20.0
    assert d["redundant_contigs"] is False and d["mito_size_estimate"] == 8080 and d["longest_contig"]["name"] == "mt"
    assert len(d["missing_core_genes"]) == 14 and d["heteroplasmy"] is None and "1/15 core genes" in r.stdout
    gff = open(tmp_path / "o.gff").read()
    assert "ID=cox1;Name=cox1;introns=0;ref=Af293" in gff
    empty = tmp_path / "e.fa"; empty.write_text("")
    r2 = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "mito_annotate.py"), "--sample", "S1", "--mito", str(empty), "--out-json", str(tmp_path / "o2.json"), "--out-gff", str(tmp_path / "o2.gff")], capture_output=True, text=True)
    assert r2.returncode == 0 and json.load(open(tmp_path / "o2.json"))["mito_present"] is False


def test_redundant_fragments_of_the_circle_are_recognised(tmp_path, helpers):
    seqs = {"f1": helpers["random_dna"](23000, seed=11), "f2": helpers["random_dna"](23000, seed=12), "f3": helpers["random_dna"](29000, seed=13)}
    fa = helpers["write_fasta"](tmp_path / "m.fa", seqs)
    lines = []
    for i, gene in enumerate(["cox1", "cob", "nad1", "nad2"]):
        for c in ("f1", "f2", "f3"):
            lines.append(hit(f"{gene}|A|Af293", c, 1000 + 2000 * i, 2000 + 2000 * i, qs=1, qe=300, bits=400, qlen=300))
    lines.append(hit("atp6|A|Af293", "f3", 15000, 15700, qs=1, qe=250, bits=300, qlen=250))
    t = tmp_path / "g.tsv"; t.write_text("\n".join(lines) + "\n")
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "mito_annotate.py"), "--sample", "S1", "--mito", fa, "--tblastn", str(t),
                        "--out-json", str(tmp_path / "o.json"), "--out-gff", str(tmp_path / "o.gff")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(tmp_path / "o.json"))
    assert d["redundant_contigs"] is True and d["mito_size"] == 75000 and d["mito_size_estimate"] == 29000 and d["longest_contig"]["name"] == "f3"
    assert d["genes_with_extra_copies"] == [] and d["genes"]["cox1"]["contigs_with_gene"] == 3 and d["genes"]["atp6"]["contigs_with_gene"] == 1
    assert "redundant" in r.stdout and "redundant fragments" in d["note"]
