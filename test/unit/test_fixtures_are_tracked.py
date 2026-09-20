"""Every test fixture on disk must be committed, or CI tests a file that only exists locally.

`.gitignore` excludes `test/reads/*` and re-admits a few filename patterns. A fixture whose name
matches none of them is skipped silently by `git add -A`: the suite passes on the machine that
created it and fails in CI with "file not found", which reads like a broken test rather than a
missing commit. That is exactly how the RNA-seq routing fixtures were lost.

Skipped when git is unavailable or this is not a work tree, so a source tarball still tests clean.
"""
import os
import subprocess

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Directories whose contents are inputs to tests, so an untracked file there is a landmine.
FIXTURE_DIRS = ["test/reads", "test/fixtures", "test/expected", "test/controls", "test/mini_db"]
IGNORE_NAMES = {".DS_Store"}


def git(*args):
    return subprocess.run(["git", "-C", REPO, *args], capture_output=True, text=True)


def is_work_tree():
    r = git("rev-parse", "--is-inside-work-tree")
    return r.returncode == 0 and r.stdout.strip() == "true"


pytestmark = pytest.mark.skipif(not is_work_tree(), reason="not a git work tree")


def tracked_in(d):
    r = git("ls-files", d)
    return {line.strip() for line in r.stdout.splitlines() if line.strip()}


def on_disk_in(d):
    out = set()
    full = os.path.join(REPO, d)
    for root, dirs, files in os.walk(full):
        dirs[:] = [x for x in dirs if x != "__pycache__"]
        for f in files:
            if f in IGNORE_NAMES or f.endswith((".pyc", ".pyo")):
                continue
            out.add(os.path.relpath(os.path.join(root, f), REPO))
    return out


@pytest.mark.parametrize("d", FIXTURE_DIRS)
def test_every_fixture_on_disk_is_tracked(d):
    if not os.path.isdir(os.path.join(REPO, d)):
        pytest.skip(f"{d} not present")
    missing = sorted(on_disk_in(d) - tracked_in(d))
    assert not missing, (
        "untracked fixture(s) under %s — the suite passes here and fails in CI: %s. "
        "Add them (git add -f if .gitignore excludes the name) and widen the .gitignore exception."
        % (d, ", ".join(missing)))
