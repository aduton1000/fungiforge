"""Unit tests for bin/control_genotype.py — read-level verification of a control's genotype.
Pure-Python SAM evaluation is tested with hand-written records; the minimap2 integration test
simulates reads from the wild-type locus and from a TR34/L98H-mutated copy and is skipped when
minimap2 is not on PATH (it is installed in CI)."""
import json, os, random, shutil, subprocess, sys
import pytest
import control_genotype as cg  # noqa: E402  (sys.path set in conftest)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SPEC = os.path.join(ROOT, "test", "controls", "afum_cyp51A_controls.json")
CHECKS = [{"name": "TR34", "type": "insertion", "window": [100, 133], "min_len": 30, "max_len": 40},
          {"name": "L98H", "type": "snv", "pos": 400, "ref": "T", "alt": "A"}]


def sam(qname, flag, pos, cigar, seq, mapq=60):
    return "\t".join([qname, str(flag), "locus", str(pos), str(mapq), cigar, "*", "0", "0", seq, "*"])


def test_cigar_walk_tracks_reference_positions():
    ops = cg.parse_cigar("3M2I2M1D2M")
    steps = list(cg.walk(10, ops, "AAAGGCCTT"))
    assert steps[:3] == [(10, "M", 1, "A"), (11, "M", 1, "A"), (12, "M", 1, "A")]
    assert steps[3] == (12, "I", 2, "GG")
    assert steps[4:6] == [(13, "M", 1, "C"), (14, "M", 1, "C")]
    assert steps[6] == (15, "D", 1, "")
    assert steps[7:] == [(16, "M", 1, "T"), (17, "M", 1, "T")]


def test_evaluate_reads_counts_insertions_and_bases():
    lines = [
        sam("r1", 0, 50, "70M34I200M", "A" * 304),             # spans window, 34-bp insertion at ref 119 -> TR34 present
        sam("r2", 0, 50, "270M", "A" * 270),                    # spans window, no insertion
        sam("r3", 0, 110, "200M", "A" * 200),                   # starts inside the window -> not spanning
        sam("r4", 0, 350, "60M", "T" * 50 + "A" + "T" * 9),     # base at 400 is A (alt)
        sam("r5", 0, 350, "60M", "T" * 60),                     # base at 400 is T (ref)
        sam("r6", 256, 350, "60M", "T" * 60),                   # secondary: ignored
        sam("r7", 0, 350, "60M", "T" * 60, mapq=0),             # mapq 0: ignored
        sam("r8", 0, 50, "70M12I200M", "A" * 282),              # insertion too short for TR34 (12 bp)
    ]
    obs = cg.evaluate_reads(lines, CHECKS)
    assert obs["TR34"]["spanning"] == 3 and obs["TR34"]["with_insertion"] == 1 and obs["TR34"]["insertion_lengths"] == {34: 1}
    assert obs["L98H"]["depth"] == 2 and obs["L98H"]["bases"] == {"A": 1, "T": 1} and obs["L98H"]["alt_frac"] == 0.5


def test_call_thresholds():
    obs = {"TR34": {"type": "insertion", "spanning": 40, "frac": 0.95},
           "TR46": {"type": "insertion", "spanning": 40, "frac": 0.0},
           "L98H": {"type": "snv", "depth": 30, "alt_frac": 0.5},
           "low":  {"type": "snv", "depth": 3, "alt_frac": 1.0}}
    assert cg.call(obs, 10, 0.8) == {"TR34": "present", "TR46": "absent", "L98H": "mixed", "low": "insufficient"}


def test_spec_file_is_consistent_with_its_locus():
    with open(SPEC) as fh:
        spec = json.load(fh)
    locus_path = os.path.join(os.path.dirname(SPEC), spec["locus"])
    seq = "".join(l.strip() for l in open(locus_path) if not l.startswith(">"))
    tr = [c for c in spec["checks"] if c["name"] == "TR34"][0]
    unit = seq[tr["window"][0] - 1:tr["window"][1]]
    assert len(unit) == 34 and seq.count(unit) == 1                 # wild-type: exactly one copy of the TR34 unit
    l98 = [c for c in spec["checks"] if c["name"] == "L98H"][0]
    assert seq[l98["pos"] - 1] == l98["ref"] == "T"
    assert seq[l98["pos"] - 2:l98["pos"] + 1] == "CTC"               # codon 98 = Leu (CTC) in the wild-type
    assert seq[1852:1855] == "ATG"                                   # canonical start codon
    assert set(spec["controls"]) >= {"TR34_L98H", "wild_type"}


def simulate_reads(seq, n, length, seed, err=0.02):
    rng = random.Random(seed); out = []
    for i in range(n):
        s = rng.randint(0, max(0, len(seq) - length)); frag = list(seq[s:s + length])
        for j in range(len(frag)):
            if rng.random() < err:
                frag[j] = rng.choice("ACGT")
        out.append(f"@r{i}\n{''.join(frag)}\n+\n{'I' * len(frag)}\n")
    return "".join(out)


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not on PATH")
def test_minimap2_integration_calls_tr34_l98h_and_wild_type(tmp_path):
    with open(SPEC) as fh:
        spec = json.load(fh)
    wt = "".join(l.strip() for l in open(os.path.join(os.path.dirname(SPEC), spec["locus"])) if not l.startswith(">"))
    tr = [c for c in spec["checks"] if c["name"] == "TR34"][0]; l98 = [c for c in spec["checks"] if c["name"] == "L98H"][0]
    unit = wt[tr["window"][0] - 1:tr["window"][1]]
    mut = wt[:tr["window"][1]] + unit + wt[tr["window"][1]:]           # duplicate the unit in tandem
    p = l98["pos"] - 1 + 34                                             # L98 shifts by one unit length
    assert mut[p] == "T"; mut = mut[:p] + "A" + mut[p + 1:]
    r_mut = tmp_path / "mut.fq"; r_mut.write_text(simulate_reads(mut, 400, 1500, 1))
    r_wt = tmp_path / "wt.fq"; r_wt.write_text(simulate_reads(wt, 400, 1500, 2))
    for reads, control, expect_pass in ((r_mut, "TR34_L98H", True), (r_wt, "wild_type", True), (r_mut, "wild_type", False)):
        out = tmp_path / f"{reads.stem}_{control}.json"
        r = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "control_genotype.py"), "--spec", SPEC, "--control", control,
                            "--reads", str(reads), "--platform", "ont", "--out", str(out)], capture_output=True, text=True)
        res = json.load(open(out))
        assert (r.returncode == 0) == expect_pass, r.stdout + r.stderr
        assert res["passed"] == expect_pass
    res = json.load(open(tmp_path / "mut_TR34_L98H.json"))
    assert res["calls"] == {"TR34": "present", "TR46": "absent", "L98H": "present"}
    lengths = {int(k): v for k, v in res["observations"]["TR34"]["insertion_lengths"].items()}
    assert 32 <= max(lengths, key=lengths.get) <= 36                      # the modal insertion is the 34-bp unit
