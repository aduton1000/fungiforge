"""Unit tests for bin/protect_organelle.py — organelle contigs survive haplotig purging.

purge_dups judges contigs by read depth against the nuclear peak. A mitochondrial genome sits far
above that peak and is dropped as collapsed duplication: on a real A. fumigatus assembly it removed
nine of ten organelle-shaped contigs including the 30.7 kb, 25.5 % GC mitogenome, which emptied
every mitochondrial column downstream. Only contigs carrying core mitochondrial genes come back.
"""
import json
import os
import subprocess
import sys

import protect_organelle as po  # noqa: E402  (sys.path set in conftest)

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(os.path.dirname(os.path.dirname(HERE)), "bin")


def fasta(path, seqs):
    with open(path, "w") as fh:
        for n, s in seqs.items():
            fh.write(">%s\n%s\n" % (n, s))
    return str(path)


def hits(path, rows):
    """tblastn outfmt 6: qseqid sseqid ... evalue bitscore."""
    with open(path, "w") as fh:
        for q, s in rows:
            fh.write("\t".join([q, s, "90", "100", "1", "100", "1", "100", "1e-30", "200"]) + "\n")
    return str(path)


MITO = "AT" * 15000          # 30 kb, 0 % GC — well under the GC ceiling
NUC = "GCAT" * 30000         # 50 % GC


def test_a_purged_mitochondrial_contig_is_restored(tmp_path):
    asm = fasta(tmp_path / "a.fa", {"contig_69": MITO, "contig_1": NUC, "contig_2": NUC})
    pur = fasta(tmp_path / "p.fa", {"contig_1_1": NUC, "contig_2_1": NUC})
    tb = hits(tmp_path / "h.tsv", [("cox1", "contig_69"), ("cob", "contig_69")])
    pre, post, mito, dropped = po.restore(asm, pur, tb)
    assert mito == ["contig_69"] and dropped == ["contig_69"]


def test_a_renamed_surviving_contig_is_not_duplicated(tmp_path):
    """get_seqs rewrites contig_69 to contig_69_1; that is the same contig, not a missing one."""
    asm = fasta(tmp_path / "a.fa", {"contig_69": MITO, "contig_1": NUC})
    pur = fasta(tmp_path / "p.fa", {"contig_69_1": MITO, "contig_1_1": NUC})
    tb = hits(tmp_path / "h.tsv", [("cox1", "contig_69"), ("cob", "contig_69")])
    _pre, _post, mito, dropped = po.restore(asm, pur, tb)
    assert mito == ["contig_69"] and dropped == []


def test_a_haplotig_without_core_genes_stays_purged(tmp_path):
    asm = fasta(tmp_path / "a.fa", {"contig_5": NUC, "contig_1": NUC})
    pur = fasta(tmp_path / "p.fa", {"contig_1_1": NUC})
    tb = hits(tmp_path / "h.tsv", [])
    _pre, _post, mito, dropped = po.restore(asm, pur, tb)
    assert mito == [] and dropped == []


def test_one_core_gene_is_not_enough(tmp_path):
    asm = fasta(tmp_path / "a.fa", {"contig_69": MITO, "contig_1": NUC})
    pur = fasta(tmp_path / "p.fa", {"contig_1_1": NUC})
    tb = hits(tmp_path / "h.tsv", [("cox1", "contig_69")])
    _pre, _post, mito, dropped = po.restore(asm, pur, tb)
    assert mito == [] and dropped == []


def test_a_long_high_gc_numt_is_not_restored(tmp_path):
    """Nuclear copies of mitochondrial genes must not drag a chromosome back in."""
    chrom = NUC * 30          # 3.6 Mb at 50 % GC
    asm = fasta(tmp_path / "a.fa", {"contig_13": chrom, "contig_1": NUC})
    pur = fasta(tmp_path / "p.fa", {"contig_1_1": NUC})
    tb = hits(tmp_path / "h.tsv", [("cox1", "contig_13"), ("cob", "contig_13")])
    _pre, _post, mito, dropped = po.restore(asm, pur, tb)
    assert mito == [] and dropped == []


def test_cli_writes_the_union_and_a_report(tmp_path):
    asm = fasta(tmp_path / "a.fa", {"contig_69": MITO, "contig_1": NUC})
    pur = fasta(tmp_path / "p.fa", {"contig_1_1": NUC})
    tb = hits(tmp_path / "h.tsv", [("cox1", "contig_69"), ("cob", "contig_69")])
    out, rep = tmp_path / "r.fa", tmp_path / "r.json"
    r = subprocess.run([sys.executable, os.path.join(BIN, "protect_organelle.py"),
                        "--assembly", asm, "--purged", pur, "--tblastn", tb,
                        "--out", str(out), "--json", str(rep)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    names = [l[1:].strip() for l in open(out) if l.startswith(">")]
    assert names == ["contig_1_1", "contig_69"]
    doc = json.load(open(rep))
    assert doc["n_restored"] == 1 and doc["restored"] == ["contig_69"] and doc["n_contigs_final"] == 2
