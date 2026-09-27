"""Static check: a module branch that runs no tool must say so with a not-applicable skip.

`ff_finalize` reports `skipped` for a stage that ran nothing and only recorded plain skips, and
the master table counts `skipped` in stages_failed (that is deliberate: a stage that wrote an
empty table because its database was absent must not read as success). A branch that by design
runs no tool — nothing to polish in an Illumina-only assembly, the passthrough is the whole job —
must record its skip with `ff_skip --not-applicable`, or every such isolate is flagged as failed.
DF-005 on the dev deployment (2026-09-21) showed `polish:skipped` for exactly this reason.
Stub tests never run real script blocks, so only a static check catches it before a real run.
"""
import glob
import os
import re

import pytest

from test_module_shell_vars import script_blocks

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MODULES = sorted(glob.glob(os.path.join(REPO, "modules", "*.nf")))

RUN = re.compile(r"\bff_run\b")
SKIP = re.compile(r"\bff_skip\s+(--not-applicable\s+)?(\S+)")
STAGE_SKIPPED = re.compile(r"--stage-skipped")


def plain_skip_only_blocks(text):
    """Branches with ff_init, no ff_run, no explicit --stage-skipped, and at least one plain ff_skip."""
    bad = []
    for block in script_blocks(text):
        if "ff_init" not in block or RUN.search(block) or STAGE_SKIPPED.search(block):
            continue
        plain = [m.group(2) for m in SKIP.finditer(block) if not m.group(1)]
        if plain:
            bad.append(plain)
    return bad


@pytest.mark.parametrize("path", MODULES, ids=[os.path.basename(p) for p in MODULES])
def test_no_branch_only_plain_skips(path):
    bad = plain_skip_only_blocks(open(path).read())
    assert not bad, (
        "%s has a branch that runs no tool and only plain-skips %s: it will report `skipped` and "
        "land in stages_failed. If the step does not apply to the isolate by design, use "
        "`ff_skip --not-applicable`; if it is a real omission, the branch is right to report it."
        % (os.path.basename(path), bad))


def test_check_detects_the_df005_shape():
    block = ("\n  script:\n    source ff_status.sh; ff_init polish S x.json\n"
             "    ff_skip polypolish \"nothing to polish\"\n    cp a b\n    ff_finalize\n  stub:\n")
    assert plain_skip_only_blocks(block) == [["polypolish"]]


def test_check_accepts_not_applicable_and_real_runs():
    ok1 = ("\n  script:\n    ff_init polish S x.json\n    ff_skip --not-applicable polypolish \"n/a\"\n"
           "    ff_finalize\n  stub:\n")
    ok2 = ("\n  script:\n    ff_init x S x.json\n    ff_skip genomad \"db absent\"\n    ff_run t -- true\n"
           "    ff_finalize\n  stub:\n")
    assert plain_skip_only_blocks(ok1) == [] and plain_skip_only_blocks(ok2) == []
