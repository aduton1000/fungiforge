#!/usr/bin/env python3
"""Emit provenance.json — tool versions, container digests, git commit — for a run
(forge reproducibility invariant). Best-effort: records whatever is resolvable."""
from __future__ import annotations
import argparse, json, os, subprocess


def sh(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="provenance.json")
    ap.add_argument("--repo", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    a = ap.parse_args()
    prov = {
        "pipeline": "fungiforge", "version": "0.1.0",
        "git_commit": sh(f"git -C {a.repo} rev-parse HEAD"),
        "git_dirty": bool(sh(f"git -C {a.repo} status --porcelain")),
        "nextflow": sh("nextflow -version 2>&1 | head -1"),
        "docker": sh("docker --version"),
        "host": sh("uname -a"),
        "tools": {t: sh(f"{t} --version 2>&1 | head -1") for t in
                  ("flye", "medaka_consensus", "funannotate", "antismash", "sourmash", "skani")},
    }
    json.dump(prov, open(a.out, "w"), indent=2)
    print(f"[provenance] -> {a.out}")


if __name__ == "__main__":
    main()
