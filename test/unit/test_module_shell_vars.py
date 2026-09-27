"""Static check: no module task script may reference an unset shell variable.

Task scripts run under `bash -ue`, so a reference to a variable nothing assigns kills the task
at that line. Stub tests cannot catch it because a stub never executes the real script block, and
a run does not reveal it until the stage is reached, which for assembly or annotation is hours in.
One such leftover (`$PURGED` in stage 02, orphaned when haplotig purging moved to stage 02c)
failed a cluster run only after Flye had finished a 30 Mb assembly.

The parse is deliberately conservative: shell references are the backslash-escaped `\\$NAME` /
`\\${NAME}` forms, since an unescaped `${...}` is Groovy interpolation resolved before bash sees it.
"""
import glob
import os
import re

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Set by the shell, the container image, or conf/base.config's environment — not by the script.
SAFE = set("""
PWD OLDPWD HOME PATH TMPDIR TMP USER LOGNAME HOSTNAME SHELL IFS RANDOM LINENO SECONDS
PS1 PS2 BASH BASH_VERSION BASHPID PPID UID EUID LANG LC_ALL TERM
FF_RC FF_STAGE FF_TOOLS FF_SAMPLE FF_JSON FF_BEST_EFFORT FF_T0
NXF_TASK_WORKDIR NXF_DEBUG NXF_OPTS
APPTAINER_BIND SINGULARITY_BIND OMP_NUM_THREADS JAVA_TOOL_OPTIONS GENEMARK_PATH
AUGUSTUS_CONFIG_PATH FUNANNOTATE_DB EGGNOG_DATA_DIR SIGNALP_DIR PERL5LIB PYTHONPATH
""".split())

ASSIGN = re.compile(
    r"""(?:^|[;&|(){}\s])(?:local\s+|export\s+|declare\s+(?:-\w+\s+)?|readonly\s+)?"""
    r"""([A-Za-z_][A-Za-z0-9_]*)=(?!=)""", re.M)
LOOP = re.compile(r"\bfor\s+([A-Za-z_][A-Za-z0-9_]*)\s+in\b"
                  r"|\bread\s+(?:-\w+\s+)*([A-Za-z_][A-Za-z0-9_]*)")
REF = re.compile(r"\\\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?")
DEFAULTED = re.compile(r"\\\$\{([A-Za-z_][A-Za-z0-9_]*)(?::-|:=|:\?|\+|-)")


def script_blocks(text):
    """The shell body of each `script:` block (a `stub:` block is not run in production)."""
    return re.findall(r"\n\s*script:\s*\n(.*?)\n\s*(?:stub:|\}\s*$)", text, re.S)


def unbound(block):
    defined = set(ASSIGN.findall(block)) | SAFE | set(DEFAULTED.findall(block))
    for m in LOOP.finditer(block):
        defined |= {g for g in m.groups() if g}
    return sorted({m.group(1) for m in REF.finditer(block) if m.group(1) not in defined})


MODULES = sorted(glob.glob(os.path.join(REPO, "modules", "*.nf")))


def test_modules_are_discovered():
    assert len(MODULES) > 20, "module glob found almost nothing — the path is wrong"


@pytest.mark.parametrize("path", MODULES, ids=[os.path.basename(p) for p in MODULES])
def test_no_unbound_shell_variables(path):
    for block in script_blocks(open(path).read()):
        names = unbound(block)
        assert not names, (
            "%s references shell variable(s) nothing assigns: %s. Under `bash -ue` the task dies "
            "at that line. Assign it, give it a default (${NAME:=...}), or add it to SAFE if the "
            "image or the config exports it." % (os.path.basename(path), ", ".join(names)))


def test_check_detects_a_planted_reference():
    """Guard the guard: the exact shape that broke a cluster run must still be caught."""
    block = 'X=1\nprintf "%s" "\\$X" "\\$PURGED" > out.json\n'
    assert unbound(block) == ["PURGED"]


def test_check_accepts_assignments_loops_and_defaults():
    block = ('N=0\nfor f in a b; do echo "\\$f"; done\n'
             'export Q="\\$PWD/tmp"\n: "\\${R:=fallback}"\necho "\\$N \\$Q \\$R"\n')
    assert unbound(block) == []
