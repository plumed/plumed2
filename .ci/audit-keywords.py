#!/usr/bin/env python3
"""Keywords audit: compare the keywords that the actions register (json/syntax.json, the source of the
manual) with the keywords that their code reads.

Two checks:
  NOT READ    a registered keyword is added to an example of the manual that uses the action
              (a flag as it is, a keyword with its default or with 1) and PLUMED answers
              "cannot understand the following words from the input line": no code reads it,
              so the manual documents a keyword that every input using it is refused with.
              Inputs that end with a crash, an assertion or an uncaught exception are reported as in
              audit-inputs.py.
  NOT REGISTERED  (static) parse("KEY") in a file that registers actions, for a KEY that none of the
              actions of that file has in json/syntax.json. Reading such a keyword stops at an assertion.
              The keywords read only after keys.exists("KEY") are listed apart: they can be correct.

Usage (from the root of a compiled source tree with json/syntax.json, after "source sourceme.sh"):
  python3 .ci/audit-keywords.py [--jobs N] [--only REGEX] [--static-only]
REGEX selects the actions by name. Each action is tried in the first example of the manual where it
appears on a single line and that runs without errors as written.
"""

import argparse
import concurrent.futures
import glob
import importlib.util
import json
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "audit_inputs", os.path.join(HERE, "audit-inputs.py")
)
ai = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ai)
ROOT = ai.ROOT

SKIP_TYPES = ("hidden", "deprecated")
NOT_UNDERSTOOD = "cannot understand the following words from the input line"


def static_check(syntax):
    parse_re = re.compile(
        r'\bparse(?:Flag|Vector|Numbered|NumberedVector|AtomList|ArgumentList)?\s*\(\s*"([A-Za-z0-9_]+)"'
    )
    exists_re = re.compile(r'exists\s*\(\s*"([A-Za-z0-9_]+)"')
    register_re = re.compile(
        r'PLUMED_REGISTER_ACTION\s*\(\s*\w+\s*,\s*"([A-Za-z0-9_]+)"'
    )
    found = []
    for fname in sorted(
        glob.glob(ROOT + "/src/*/*.cpp") + glob.glob(ROOT + "/src/*/*.h")
    ):
        if fname.startswith((ROOT + "/src/include/", ROOT + "/src/lib/")):
            continue
        with open(fname, errors="replace") as f:
            text = f.read()
        # comments are not code (the manual examples are in comments)
        text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
        text = re.sub(r"//[^\n]*", "", text)
        actions = [a for a in register_re.findall(text) if a in syntax]
        if not actions:
            continue
        known = set()
        for a in actions:
            known.update(syntax[a]["syntax"].keys())
        guarded = set(exists_re.findall(text))
        for key in sorted(set(parse_re.findall(text))):
            if key in known:
                continue
            kind = (
                "NOT REGISTERED (after exists)" if key in guarded else "NOT REGISTERED"
            )
            found.append((kind, os.path.relpath(fname, ROOT), key, ",".join(actions)))
    return found


def keyword_token(key, entry):
    if entry["type"] == "flag":
        return key
    default = entry.get("default", "")
    if default and default != "off" and " " not in default:
        return f"{key}={default}"
    return f"{key}=1"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--jobs", type=int, default=os.cpu_count())
    parser.add_argument(
        "--only", default="", help="regular expression on the action names"
    )
    parser.add_argument("--static-only", action="store_true")
    parser.add_argument("--natoms", type=int, default=100000)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--plumed", default="plumed")
    args = parser.parse_args()
    args.run = False

    with open(os.path.join(ROOT, "json", "syntax.json")) as f:
        syntax = {
            k: v
            for k, v in json.load(f).items()
            if isinstance(v, dict) and "syntax" in v
        }

    nfound = 0
    for kind, fname, key, actions in static_check(syntax):
        if re.search(args.only, actions):
            print(f"{kind:30s} {fname}: {key} (actions {actions})")
            nfound += 1
    if args.static_only:
        return 1 if nfound else 0

    # candidate examples: single-line uses of each action
    candidates = {}
    for name, settings, lines in ai.extract_examples():
        if ai.skip_reason(settings, lines):
            continue
        for i, line in enumerate(lines):
            code = line.split("#")[0]
            if "..." in code or not code.strip():
                continue
            m = ai.ACTION.match(code)
            if (
                not m
                or m.group(1) not in syntax
                or not re.search(args.only, m.group(1))
            ):
                continue
            candidates.setdefault(m.group(1), []).append((name, settings, lines, i))

    args.data = tempfile.mkdtemp(prefix="auditdata")
    shutil.copytree(
        os.path.join(ROOT, "regtest"),
        os.path.join(args.data, "regtest"),
        symlinks=True,
        ignore=shutil.ignore_patterns("tmp"),
    )
    shutil.copytree(
        os.path.join(ROOT, "new-manual", "extras"), os.path.join(args.data, "extras")
    )

    def works(c):
        kind, output = ai.run_one(c[2], c[1], args)
        return kind is None and "ERROR" not in output

    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        # the first example (of up to three) that runs as written
        chosen = {}
        for action, cands in candidates.items():
            for c in cands[:3]:
                if pool.submit(works, c).result():
                    chosen[action] = c
                    break
        untested = sorted(set(syntax) - set(chosen))
        runs = []
        for action, (name, settings, lines, i) in sorted(chosen.items()):
            code = lines[i].split("#")[0]
            for key, entry in syntax[action]["syntax"].items():
                if not isinstance(entry, dict) or "type" not in entry:
                    continue
                if entry["type"] in SKIP_TYPES or key == "LABEL":
                    continue
                if re.search(r"\b" + key + r"[0-9]*(=|\s|$)", code):
                    continue
                new = (
                    lines[:i]
                    + [code.rstrip() + " " + keyword_token(key, entry)]
                    + lines[i + 1 :]
                )
                runs.append((action, key, name, settings, new))
        print(
            f"{len(chosen)} actions tried in an example, {len(untested)} without a working single-line example, {len(runs)} runs",
            flush=True,
        )
        futures = {pool.submit(ai.run_one, r[4], r[3], args): r for r in runs}
        results = []
        for fut in concurrent.futures.as_completed(futures):
            action, key, name, settings, lines = futures[fut]
            kind, output = fut.result()
            if kind is None:
                for line in output.splitlines():
                    if NOT_UNDERSTOOD in line and re.search(
                        r"\b" + key + r"\b", line.split(NOT_UNDERSTOOD)[1]
                    ):
                        kind = "NOT READ"
                        break
            if kind:
                message = ""
                for line in output.splitlines():
                    if NOT_UNDERSTOOD in line or any(
                        p.search(line) for _, p in ai.SIGNATURES
                    ):
                        message = line.strip()[:160]
                        break
                results.append((kind, action, key, name, message))
    shutil.rmtree(args.data, ignore_errors=True)
    for kind, action, key, name, message in sorted(results):
        print(f"{kind:30s} {action} {key} ({name}): {message}")
    nfound += len(results)
    print(f"{nfound} findings")
    return 1 if nfound else 0


if __name__ == "__main__":
    sys.exit(main())
