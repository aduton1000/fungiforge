"""Shared fixtures for the fungiforge unit tests (synthetic sequences, GenBank writer)."""
import os, random, sys
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
BIN = os.path.join(ROOT, "bin")
for p in (ROOT, BIN):
    if p not in sys.path:
        sys.path.insert(0, p)


def random_dna(n, seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def random_protein(n, seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n))


# the canonical A. fumigatus TR34 unit is 34 bp; any 34-mer duplicated in tandem is what the
# structural scanner is meant to find, so a synthetic unit is a faithful test
TR34_UNIT = "GAATCACGCGGTCCGGATGTGTGCTGAGCCGAAT"
assert len(TR34_UNIT) == 34


def write_fasta(path, seqs, width=0):
    with open(path, "w") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n")
            if width:
                for i in range(0, len(s), width):
                    fh.write(s[i:i + width] + "\n")
            else:
                fh.write(s + "\n")
    return str(path)


def write_gbk(path, seq, features, record_id="contig_1"):
    """features: list of (start0, end0, strand, type, qualifiers dict). Needs Biopython."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    rec = SeqRecord(Seq(seq), id=record_id, name=record_id, description="synthetic",
                    annotations={"molecule_type": "DNA"})
    for s, e, strand, ftype, quals in features:
        rec.features.append(SeqFeature(FeatureLocation(s, e, strand=strand), type=ftype, qualifiers=quals))
    SeqIO.write(rec, str(path), "genbank")
    return str(path)


@pytest.fixture
def helpers():
    return {"random_dna": random_dna, "random_protein": random_protein, "write_fasta": write_fasta,
            "write_gbk": write_gbk, "TR34_UNIT": TR34_UNIT, "ROOT": ROOT, "BIN": BIN}
