"""Unit tests for bin/read_genotype.py — CDS codon coordinates from a GenBank (both strands,
spliced), per-read codon translation, zygosity, indel support at a window, and the CLI."""
import json, os, subprocess, sys
import pytest

pytest.importorskip("Bio")
import read_genotype as rg  # noqa: E402  (sys.path set in conftest)

CONTIG = "contig_1"
FIRST = {aa: c for c, aa in sorted(rg.CODON.items(), reverse=True)}   # one codon per amino acid


def encode(protein):
    return "".join(FIRST[a] for a in protein)


def write_gbk(path, seq, cds_parts, strand, protein, locus="FUN_000010"):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    rec = SeqRecord(Seq(seq), id=CONTIG, name=CONTIG, description="synthetic", annotations={"molecule_type": "DNA"})
    locs = [FeatureLocation(s, e, strand=strand) for s, e in cds_parts]
    loc = CompoundLocation(locs) if len(locs) > 1 else locs[0]
    rec.features.append(SeqFeature(loc, type="CDS", qualifiers={"locus_tag": [locus], "protein_id": [f"ncbi:{locus}-T1"], "translation": [protein], "codon_start": ["1"]}))
    SeqIO.write(rec, str(path), "genbank")
    return str(path)


def sam_read(name, pos, seq, cigar=None, flag=0, mapq=60, contig=CONTIG):
    return f"{name}\t{flag}\t{contig}\t{pos}\t{mapq}\t{cigar or f'{len(seq)}M'}\t*\t0\t0\t{seq}\t*"


def build_gene(helpers, strand=1):
    """protein of 40 aa; exon1 = codons 1-20 at 1001.., intron 60 bp, exon2 = codons 21-40."""
    prot = helpers["random_protein"](40, seed=9).replace("*", "A")
    cds = encode(prot)
    ex1, ex2 = cds[:60], cds[60:]
    left, intron, right = helpers["random_dna"](1000, seed=1), "GT" + helpers["random_dna"](56, seed=2) + "AG", helpers["random_dna"](1000, seed=3)
    if strand == 1:
        seq = left + ex1 + intron + ex2 + right
        parts = [(1000, 1060), (1120, 1180)]
    else:
        seq = left + rg.revcomp(ex2) + rg.revcomp(intron) + rg.revcomp(ex1) + right     # gene on the minus strand
        parts = [(1000, 1060), (1120, 1180)]
    return prot, seq, parts


def test_cds_codon_coords_plus_and_minus_strand_spliced(tmp_path, helpers):
    for strand in (1, -1):
        prot, seq, parts = build_gene(helpers, strand)
        gbk = write_gbk(tmp_path / f"g{strand}.gbk", seq, parts, strand, prot)
        cds = rg.cds_codon_coords(gbk, "FUN_000010-T1", prot)
        assert cds["contig"] == CONTIG and cds["strand"] == strand and cds["n_residues"] == 40
        # every codon, read from the genome in transcript orientation, translates to the protein
        got = ""
        for trip in cds["codons"]:
            codon = "".join(seq[p - 1] for p in trip)
            got += rg.translate_codon(codon if strand == 1 else rg.complement(codon))   # trip is in transcript order
        assert got == prot
        if strand == 1:
            assert cds["codons"][19][-1] == 1060 and cds["codons"][20][0] == 1121      # exon junction
        else:
            assert cds["codons"][0] == (1180, 1179, 1178) and cds["codons"][20][0] == 1060
    # match by name when the translation is not given; None for an unknown protein
    assert rg.cds_codon_coords(gbk, "FUN_000010-T1")["n_residues"] == 40
    assert rg.cds_codon_coords(gbk, "FUN_999999-T1") is None


def reads_over(seq, start, end, n, step=7, length=120):
    out = []
    for i in range(n):
        p = start - 60 + (i * step) % 40
        out.append(sam_read(f"r{i}", p, seq[p - 1:p - 1 + length]))
    return out


def test_codon_counts_zygosity_and_junction_codon(tmp_path, helpers):
    prot, seq, parts = build_gene(helpers, 1)
    gbk = write_gbk(tmp_path / "g.gbk", seq, parts, 1, prot)
    cds = rg.cds_codon_coords(gbk, "FUN_000010-T1", prot)
    r10 = cds["codons"][9]                  # residue 10, exon 1
    r21 = cds["codons"][20]                 # residue 21, first codon of exon 2
    alt = FIRST["H"] if prot[9] != "H" else FIRST["K"]
    reads = []
    for i in range(30):                     # 30 reads over the region; half carry the variant at residue 10
        p = 980 + (i * 3) % 30
        s = list(seq[p - 1:p - 1 + 150])          # reads end before residue 40
        if i % 2 == 0:
            for k, base in zip(r10, alt):
                s[k - p] = base
        reads.append(sam_read(f"r{i}", p, "".join(s)))
    geno = rg.genotype_residues(("sam", reads), cds, [10, 21, 40], min_reads=10)
    assert geno[10]["call"] == "heterozygous" and set(geno[10]["alleles"]) == {prot[9], "H" if alt == FIRST["H"] else "K"} and geno[10]["depth"] == 30
    assert geno[21]["call"] == "homozygous" and geno[21]["major"] == prot[20] and geno[21]["codon"] == list(r21)
    assert geno[40]["call"] == "insufficient"                           # reads end before the last codon
    # a read with a deletion over the codon counts as 'del'; secondary/unmapped/low-MAPQ reads are ignored
    p = r10[0] - 30
    dele = sam_read("d", p, seq[p - 1:p - 1 + 30] + seq[p + 32:p + 62], cigar="30M3D30M")
    bad = [sam_read("s", p, seq[p - 1:p - 1 + 60], flag=256), sam_read("u", p, seq[p - 1:p - 1 + 60], flag=4), sam_read("q", p, seq[p - 1:p - 1 + 60], mapq=0)]
    counts = rg.codon_counts([dele] + bad, CONTIG, 1, {10: r10})
    assert counts[10] == {"del": 1}


def test_zygosity_thresholds():
    import collections
    assert rg.zygosity(collections.Counter({"L": 5}), min_reads=10)["call"] == "insufficient"
    assert rg.zygosity(collections.Counter({"L": 9, "H": 1}))["call"] == "homozygous"
    z = rg.zygosity(collections.Counter({"L": 6, "H": 4}))
    assert z["call"] == "heterozygous" and z["major"] == "L" and z["alleles"] == {"L": 0.6, "H": 0.4}
    assert rg.zygosity(collections.Counter({"L": 7, "H": 1, "X": 1, "del": 1}))["call"] == "mixed"


def test_indel_support_counts_insertions_and_deletions_in_window(helpers):
    seq = helpers["random_dna"](600, seed=5)
    win = (300, 333)
    ins = sam_read("i", 200, seq[199:332] + "A" * 34 + seq[332:460], cigar="133M34I128M")
    dele = sam_read("d", 200, seq[199:299] + seq[333:460], cigar="100M34D127M")
    plain = sam_read("p", 200, seq[199:460])
    short = sam_read("s", 290, seq[289:400])                             # does not span the window with margin
    far_ins = sam_read("f", 200, seq[199:219] + "T" * 34 + seq[219:460], cigar="20M34I241M")   # insertion outside the window
    o = rg.indel_support([ins, dele, plain, short, far_ins], CONTIG, win, 30, 40)
    assert o["spanning"] == 4 and o["with_insertion"] == 1 and o["with_deletion"] == 1
    assert o["insertion_lengths"] == {34: 1} and o["deletion_lengths"] == {34: 1} and o["ins_frac"] == 0.25 and o["del_frac"] == 0.25
    assert rg.insertion_support([ins, plain], CONTIG, win, 30, 40) == {"spanning": 2, "with_insertion": 1, "insertion_lengths": {34: 1}, "frac": 0.5}


def test_cli_with_sam(tmp_path, helpers):
    prot, seq, parts = build_gene(helpers, 1)
    gbk = write_gbk(tmp_path / "g.gbk", seq, parts, 1, prot)
    faa = helpers["write_fasta"](tmp_path / "p.faa", {"FUN_000010-T1": prot})
    sam = tmp_path / "r.sam"
    sam.write_text("@HD\tVN:1.6\n" + "\n".join(sam_read(f"r{i}", 990 + i, seq[989 + i:989 + i + 150]) for i in range(12)) + "\n")
    out = tmp_path / "g.json"
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "read_genotype.py"), "--sam", str(sam), "--gbk", gbk, "--protein", "FUN_000010-T1",
                        "--protein-fasta", faa, "--residues", "5,40", "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    doc = json.load(open(out))
    assert doc["residues"]["5"]["call"] == "homozygous" and doc["residues"]["5"]["major"] == prot[4] and doc["residues"]["40"]["call"] == "insufficient"
    assert "FUN_000010-T1\t5\thomozygous" in r.stdout


def test_gbk_to_fasta_writes_record_sequences(tmp_path, helpers):
    prot, seq, parts = build_gene(helpers, 1)
    gbk = write_gbk(tmp_path / "g.gbk", seq, parts, 1, prot)
    out = tmp_path / "records.fa"
    r = subprocess.run([sys.executable, os.path.join(helpers["BIN"], "gbk_to_fasta.py"), "--gbk", gbk, "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 0 and "1 records" in r.stdout
    lines = out.read_text().splitlines()
    assert lines[0] == f">{CONTIG}" and "".join(lines[1:]) == seq and max(len(l) for l in lines[1:]) == 80
