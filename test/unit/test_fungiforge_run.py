"""The shared launcher (share/bin/fungiforge-run) against a stub `nextflow` that prints its arguments.

Regression (2026-09-25): the launcher re-sourced the site env file on every start, and that file
resets FUNGIFORGE_WORK, so `FUNGIFORGE_WORK=~/x fungiforge-run ... -resume` in a shell that had
already loaded the env resumed nothing (a day of compute restarted from read QC). The env file
now marks itself loaded and the launcher sources it only when the shell has not.
"""
import os
import stat
import subprocess
import textwrap

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LAUNCHER = os.path.join(REPO, "share", "bin", "fungiforge-run")


def make_site(tmp_path):
    """A fake install: <root>/repo/main.nf + launcher copy, <root>/fungiforge-env.sh from the example, a stub nextflow."""
    root = tmp_path / "site"
    (root / "repo" / "share" / "bin").mkdir(parents=True)
    (root / "repo" / "main.nf").write_text("// stub\n")
    launcher = root / "repo" / "share" / "bin" / "fungiforge-run"
    launcher.write_text(open(LAUNCHER).read())
    launcher.chmod(launcher.stat().st_mode | stat.S_IEXEC)
    env = open(os.path.join(REPO, "share", "fungiforge-env.sh.example")).read()
    env = env.replace('export FUNGIFORGE_ROOT="/hpc/opt/fungiforge"', f'export FUNGIFORGE_ROOT="{root}"')
    env = env.replace('export FUNGIFORGE_DB="/hpc/data/fungiforge"', f'export FUNGIFORGE_DB="{tmp_path / "db"}"')
    assert str(root) in env and str(tmp_path / "db") in env, "example env layout changed; update the test"
    (root / "fungiforge-env.sh").write_text(env)
    (tmp_path / "db").mkdir()
    binstub = tmp_path / "stub"; binstub.mkdir()
    nf = binstub / "nextflow"
    nf.write_text("#!/bin/bash\nprintf '%s\\n' \"$@\"\n"); nf.chmod(0o755)
    return root, str(binstub)


def launch(tmp_path, root, binstub, preload_env, extra_env=None, args=("--samplesheet", "s.csv")):
    """Run the launcher from a run dir, optionally in a shell that has ALREADY sourced the env file."""
    rundir = tmp_path / "run"; rundir.mkdir(exist_ok=True)
    pre = f'source "{root}/fungiforge-env.sh"\n' if preload_env else ""
    exports = "".join(f'export {k}="{v}"\n' for k, v in (extra_env or {}).items())
    script = textwrap.dedent(f"""\
        set -e
        export PATH="{binstub}:$PATH"
        export HOME="{tmp_path}"
        {pre}{exports}cd "{rundir}"
        "{root}/repo/share/bin/fungiforge-run" {" ".join(args)}
        """)
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env={"PATH": os.environ["PATH"], "HOME": str(tmp_path)})
    assert r.returncode == 0, r.stderr
    return r.stdout.strip().split("\n"), r.stderr


def work_dir_of(argv):
    return argv[argv.index("-work-dir") + 1]


def test_env_file_marks_itself_loaded(tmp_path):
    root, _ = make_site(tmp_path)
    r = subprocess.run(["bash", "-c", f'source "{root}/fungiforge-env.sh"; echo "$FUNGIFORGE_ENV_LOADED"'], capture_output=True, text=True)
    assert os.path.realpath(r.stdout.strip()) == os.path.realpath(str(root / "fungiforge-env.sh"))


def test_explicit_work_dir_survives_when_the_shell_already_loaded_the_env(tmp_path):
    root, stub = make_site(tmp_path)
    argv, err = launch(tmp_path, root, stub, preload_env=True, extra_env={"FUNGIFORGE_WORK": str(tmp_path / "scratch")})
    assert work_dir_of(argv) == str(tmp_path / "scratch"), argv
    assert "already loaded" in err


def test_env_not_yet_loaded_is_sourced_and_clears_a_leaked_global_work_dir(tmp_path):
    """A login shell carrying another install's global FUNGIFORGE_WORK: the launcher sources the env
    file and its reset applies, so the run uses work/ in the run directory."""
    root, stub = make_site(tmp_path)
    argv, err = launch(tmp_path, root, stub, preload_env=False, extra_env={"FUNGIFORGE_WORK": str(tmp_path / "leaked")})
    assert work_dir_of(argv) == str(tmp_path / "run" / "work"), argv
    assert "already loaded" not in err


def test_nextflow_work_dir_argument_always_wins(tmp_path):
    root, stub = make_site(tmp_path)
    argv, _ = launch(tmp_path, root, stub, preload_env=True, args=("--samplesheet", "s.csv", "-work-dir", str(tmp_path / "w")))
    assert argv.count("-work-dir") == 1 and work_dir_of(argv) == str(tmp_path / "w")


def test_default_work_dir_is_the_run_directory(tmp_path):
    root, stub = make_site(tmp_path)
    argv, _ = launch(tmp_path, root, stub, preload_env=True)
    assert work_dir_of(argv) == str(tmp_path / "run" / "work") and "run" in argv and argv[-2:] == ["--samplesheet", "s.csv"]
