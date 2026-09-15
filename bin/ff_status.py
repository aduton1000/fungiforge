#!/usr/bin/env python3
"""fungiforge — stage status contract (Python side).

Every stage's JSON carries:
  "status": ok | partial | failed | skipped
  "tools":  {label: {"exit": int, "optional": bool, "seconds": int}}
  "skipped_tools": {label: reason}
  "note":   free text (set on failure/partial)

  ok       every tool that ran exited 0
  partial  an --optional tool failed and the module used its documented fallback
  failed   a required tool failed (the task exits non-zero unless the stage is best-effort)
  skipped  the stage was deliberately not run (--stage-skipped REASON)

Sub-commands
  finalize  --stage S --sample X --json F --tools T.tsv [--skips S.tsv] [--best-effort]
            [--stage-skipped REASON]   -> merges status into F (creates F if the stage's own
            writer never ran), prints one line, exits 1 if failed and not best-effort
  status    F.json [F2.json ...]        -> prints "stage<TAB>status" per file (for tests/CI)
"""
from __future__ import annotations
import argparse, json, os, sys


def read_tools(path):
    tools = {}
    if path and os.path.exists(path):
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0]:
                tools[f[0]] = {"exit": int(f[1]), "optional": f[2] == "1", "seconds": int(f[3])}
    return tools


def read_skips(path):
    skips = {}
    if path and os.path.exists(path):
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0]:
                skips[f[0]] = f[1]
    return skips


def compute_status(tools, stage_skipped=None):
    if stage_skipped:
        return "skipped", f"stage skipped: {stage_skipped}"
    failed = [k for k, v in tools.items() if v["exit"] != 0 and not v["optional"]]
    degraded = [k for k, v in tools.items() if v["exit"] != 0 and v["optional"]]
    if failed:
        return "failed", "required tool(s) failed: " + ", ".join(f"{k} (exit {tools[k]['exit']})" for k in failed)
    if degraded:
        return "partial", "optional tool(s) failed, fallback used: " + ", ".join(f"{k} (exit {tools[k]['exit']})" for k in degraded)
    return "ok", ""


def cmd_finalize(a):
    tools, skips = read_tools(a.tools), read_skips(a.skips)
    status, note = compute_status(tools, a.stage_skipped)
    doc = {}
    if os.path.exists(a.json):
        try:
            doc = json.load(open(a.json))
        except Exception:
            doc = {}
            note = (note + "; " if note else "") + "stage JSON was unreadable and has been replaced"
            if status == "ok":
                status = "failed"
    doc.setdefault("sample", a.sample)
    doc.setdefault("stage", a.stage)
    doc["status"] = status
    doc["tools"] = tools
    doc["skipped_tools"] = skips
    if note:
        doc["note"] = (doc.get("note") + " | " if doc.get("note") else "") + note
    json.dump(doc, open(a.json, "w"), indent=2)
    print(f"[ff:{a.stage}] {a.sample}: {status}" + (f" — {note}" if note else ""))
    if status == "failed" and not a.best_effort:
        sys.exit(1)


def cmd_status(a):
    for p in a.files:
        try:
            d = json.load(open(p))
            print(f"{d.get('stage', '?')}\t{d.get('status', 'unknown')}")
        except Exception:
            print(f"{os.path.basename(p)}\tunreadable")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("finalize")
    f.add_argument("--stage", required=True); f.add_argument("--sample", required=True); f.add_argument("--json", required=True)
    f.add_argument("--tools", required=True); f.add_argument("--skips"); f.add_argument("--best-effort", action="store_true")
    f.add_argument("--stage-skipped", default=None); f.set_defaults(fn=cmd_finalize)
    s = sub.add_parser("status"); s.add_argument("files", nargs="+"); s.set_defaults(fn=cmd_status)
    a = ap.parse_args(); a.fn(a)


if __name__ == "__main__":
    main()
