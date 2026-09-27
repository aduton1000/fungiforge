"""Static check: a module may not READ its stage status JSON before anything creates it.

`ff_init <stage> <sample> <file>` names the status document, but nothing writes it until
`ff_finalize` runs. A run-level stage that folds its status into a published artefact has to do
that after finalizing. Stage 16 read it before, so the cohort stage died with
`FileNotFoundError: cohort.json` after skani, MAFFT, IQ-TREE, minimap2, DIAMOND and BiG-SCAPE had
all succeeded, five hours into a run. Stub blocks never execute, so nf-test cannot catch it.

Several stages legitimately write that JSON themselves (it doubles as the stage payload), read it
back to merge a helper's output, and let ff_finalize add the status afterwards. So a read is only
flagged when nothing in the script created the file before it.
"""
import glob
import os
import re

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MODULES = sorted(glob.glob(os.path.join(REPO, "modules", "*.nf")))

READ = r"""open\(\s*["']{0}["']\s*\)|(?<![\w/])cat\s+{0}\b"""
WRITE = r"""open\(\s*["']{0}["']\s*,\s*["']w|--out-json\s+{0}\b|>\s*{0}\b|touch\s+[^\n]*{0}\b"""


def script_blocks(text):
    blocks = re.findall(r"\n\s*script:\s*\n(.*?)\n\s*(?:stub:|\}\s*$)", text, re.S)
    # Drop comment-only lines: a comment mentioning ff_finalize or the JSON name would otherwise
    # be read as the command itself, and offsets would be meaningless.
    return ["\n".join(l for l in b.splitlines() if not l.lstrip().startswith("#")) for b in blocks]


def status_file(block):
    m = re.search(r"ff_init\s+\S+\s+\S+\s+(\S+)", block)
    return m.group(1) if m else None


def first_write(block, name):
    m = re.search(WRITE.format(re.escape(name)), block)
    return m.start() if m else None


def premature_reads(block):
    name = status_file(block)
    if not name:
        return []
    # No ff_finalize at all means the status file is never written, so any read of it is wrong.
    m = re.search(r"^\s*ff_finalize\b", block, re.M)
    # No ff_finalize at all means the status file is never written, so any read of it is wrong.
    fin = m.start() if m else len(block)
    made = first_write(block, name)
    limit = fin if made is None else min(fin, made)
    return [m.start() for m in re.finditer(READ.format(re.escape(name)), block) if m.start() < limit]


def test_modules_are_discovered():
    assert len(MODULES) > 20


@pytest.mark.parametrize("path", MODULES, ids=[os.path.basename(p) for p in MODULES])
def test_status_json_is_not_read_before_it_exists(path):
    for block in script_blocks(open(path).read()):
        assert not premature_reads(block), (
            "%s reads its stage status file (%s) before anything creates it; move the block after "
            "ff_finalize." % (os.path.basename(path), status_file(block)))


def test_check_detects_a_planted_early_read():
    """Guard the guard: the exact shape that failed the cohort stage must still be caught."""
    block = ('source ff_status.sh; ff_init cohort run cohort.json --best-effort\n'
             'python3 -c \'import json; json.load(open("cohort.json"))\'\n'
             'ff_finalize\n')
    assert premature_reads(block)


def test_a_read_after_finalize_is_accepted():
    block = ('source ff_status.sh; ff_init cohort run cohort.json --best-effort\n'
             'ff_finalize\n'
             'python3 -c \'import json; json.load(open("cohort.json"))\'\n')
    assert premature_reads(block) == []


def test_write_then_read_before_finalize_is_accepted():
    """The common pattern: write the payload, read it back to merge, let ff_finalize add status."""
    block = ('source ff_status.sh; ff_init decontam s s.decontam.json\n'
             'classify --out-json s.decontam.json\n'
             'python3 -c \'import json; d = json.load(open("s.decontam.json"))\'\n'
             'ff_finalize\n')
    assert premature_reads(block) == []
