#!/usr/bin/env python3
"""Input-validation audit: run the examples of the manual, and broken variants of them, and report
the inputs that crash PLUMED instead of stopping with an error message.

Every ```plumed block of src/**/*.{cpp,h} and new-manual/**/*.md is an input. Each input is run as
it is, and then once for every variant obtained by changing one keyword of one line:
  - the keyword removed (compulsory keywords must give an error);
  - a number replaced by 0 and by -1;
  - a list of atoms or arguments replaced by its first element.
Each change of a keyword of an action is tried in the first example that has it (--all-variants: in all).
Each run is "plumed driver --parse-only" (--run adds a short trajectory and calculates).

A run is reported when it ends in a way a user should never see:
  CRASH      signal (segmentation fault, abort, FPE...) or a sanitizer report
  ASSERT     plumed_assert / plumed_massert ("+++ assertion failed")
  EXCEPTION  a std exception that is not a PLUMED error (map::at, vector::_M_range_check, stoi...)
  TIMEOUT    the run did not finish
Errors raised with error() or plumed_error() are the expected outcome and are not reported.

Usage (from the root of a compiled source tree, after "source sourceme.sh"):
  python3 .ci/audit-inputs.py [--run] [--jobs N] [--only REGEX] [--no-variants] [--out DIR]
The summary goes to stdout, the inputs and outputs of the reported runs to DIR (default audit/).
"""

import argparse
import concurrent.futures
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def extract_examples():
    """Return a list of (name, settings, lines) for every ```plumed block."""
    # src/include and src/lib only hold links to the files of the modules
    files = [f for f in sorted(glob.glob(ROOT + "/src/*/*.cpp") + glob.glob(ROOT + "/src/*/*.h"))
             if not f.startswith((ROOT + "/src/include/", ROOT + "/src/lib/"))]
    files += sorted(glob.glob(ROOT + "/new-manual/**/*.md", recursive=True))
    files += sorted(glob.glob(ROOT + "/plugins/*/README.md"))
    examples = []
    for fname in files:
        with open(fname, errors="replace") as f:
            text = f.read()
        if "```plumed" not in text:
            continue
        rel = os.path.relpath(fname, ROOT)
        for n, block in enumerate(re.findall(r"```plumed[^\n]*\n(.*?)```", text, re.S)):
            settings = {}
            lines = []
            for line in block.splitlines():
                if line.lstrip().startswith("#SETTINGS"):
                    for w in line.split()[1:]:
                        if "=" in w:
                            k, v = w.split("=", 1)
                            settings[k] = v
                    continue
                lines.append(line)
            examples.append((f"{rel}#{n}", settings, lines))
    return examples


def skip_reason(settings, lines):
    text = "\n".join(lines)
    if "__FILL__" in text:
        return "exercise with blanks"
    if int(settings.get("NREPLICAS", "1")) > 1:
        return "needs replicas"
    if re.search(r"^\s*LOAD\b", text, re.M):
        return "LOAD"
    return None


KEYVAL = re.compile(r"(\b[A-Za-z_][A-Za-z0-9_]*)=(\{[^}]*\}|[^\s]+)")
NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


ACTION = re.compile(r"^\s*(?:[^\s:]+:\s*)?([A-Za-z_][A-Za-z0-9_]*)")


def variants(lines):
    """Yield (action, description, lines) for every one-keyword change."""
    block = None
    for i, line in enumerate(lines):
        code = line.split("#")[0]
        if "..." in code:
            # first or last line of a multi-line action: the keywords in between belong to it
            if re.match(r"^\s*\.\.\.", code):
                block = None
            else:
                m = ACTION.match(code)
                block = m.group(1) if m else None
            continue
        if not code.strip():
            continue
        m = ACTION.match(code)
        action = block if block else (m.group(1) if m else "")
        for m in KEYVAL.finditer(code):
            key, val = m.group(1), m.group(2)
            if key == "LABEL":
                continue

            def changed(newtoken, what):
                new = code[: m.start()] + newtoken + code[m.end():]
                return (action, f"line {i + 1}: {what}", lines[:i] + [new] + lines[i + 1:])

            yield changed("", f"{key} removed")
            if NUMBER.match(val):
                yield changed(f"{key}=0", f"{key}=0")
                yield changed(f"{key}=-1", f"{key}=-1")
            elif "," in val and not val.startswith("{"):
                yield changed(f"{key}={val.split(',')[0]}", f"{key}={val.split(',')[0]}")


SIGNATURES = [
    ("CRASH", re.compile(r"AddressSanitizer|UndefinedBehaviorSanitizer|runtime error:|Segmentation fault|"
                         r"Floating point exception")),
    ("ASSERT", re.compile(r"\+\+\+ assertion failed")),
    ("EXCEPTION", re.compile(r"map::at|_M_range_check|std::out_of_range|std::invalid_argument|\bstoi\b|\bstod\b|"
                             r"std::bad_alloc|std::length_error|vector::reserve|basic_string::|"
                             r"terminate called|filesystem error")),
]


def classify(rc, output):
    for kind, pattern in SIGNATURES:
        if pattern.search(output):
            return kind
    if rc < 0 or rc >= 128:
        return "CRASH"
    return None


def make_trajectory(directory, natoms):
    """One frame with natoms atoms on a jittered cubic lattice (no overlapping atoms)."""
    import random
    rnd = random.Random(1)
    side = int(round(natoms ** (1.0 / 3) + 0.5))
    box = side * 0.3
    with open(os.path.join(directory, "traj.xyz"), "w") as f:
        for _ in range(2):
            f.write(f"{natoms}\n{box} {box} {box}\n")
            for n in range(natoms):
                x, y, z = n % side, (n // side) % side, n // side // side
                f.write("X %.4f %.4f %.4f\n" % tuple(0.3 * c + 0.05 * rnd.random() for c in (x, y, z)))


def run_one(lines, settings, args):
    d = tempfile.mkdtemp(prefix="audit")
    try:
        # links to the copy of the data, since some inputs write to the files they mention
        for link in ("regtest", "extras"):
            os.symlink(os.path.join(args.data, link), os.path.join(d, link))
        for f in settings.get("INPUTFILES", "").split(","):
            if f and os.path.exists(os.path.join(args.data, f)):
                dest = os.path.join(d, os.path.basename(f))
                if not os.path.exists(dest):
                    shutil.copy(os.path.join(args.data, f), dest)
        with open(os.path.join(d, "plumed.dat"), "w") as f:
            f.write("\n".join(lines) + "\n")
        cmd = [args.plumed, "--no-mpi", "driver", "--plumed", "plumed.dat", "--kt", "2.49"]
        if args.run:
            shutil.copy(args.trajectory, os.path.join(d, "traj.xyz"))
            cmd += ["--ixyz", "traj.xyz"]
        else:
            cmd += ["--parse-only", "--natoms", str(args.natoms)]
        if "MOLFILE" in settings and os.path.exists(os.path.join(args.data, settings["MOLFILE"])):
            cmd += ["--pdb", os.path.join(args.data, settings["MOLFILE"])]
        env = dict(os.environ, PLUMED_STACK_TRACE="no", PLUMED_MAXBACKUP="0")
        try:
            p = subprocess.run(cmd, cwd=d, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               timeout=args.timeout)
            output = p.stdout.decode(errors="replace")
            return classify(p.returncode, output), output
        except subprocess.TimeoutExpired as e:
            return "TIMEOUT", (e.stdout or b"").decode(errors="replace")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", action="store_true", help="calculate on a short trajectory, not only parse")
    parser.add_argument("--jobs", type=int, default=os.cpu_count())
    parser.add_argument("--only", default="", help="regular expression on the example names")
    parser.add_argument("--no-variants", action="store_true", help="run the examples as they are, only")
    parser.add_argument("--all-variants", action="store_true",
                        help="try a change in every example where it appears, not only in the first one")
    parser.add_argument("--natoms", type=int, default=100000)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--plumed", default="plumed")
    parser.add_argument("--out", default="audit")
    args = parser.parse_args()

    examples = [e for e in extract_examples() if re.search(args.only, e[0])]
    runs = []
    seen = set()
    skipped = 0
    for name, settings, lines in examples:
        if skip_reason(settings, lines):
            skipped += 1
            continue
        runs.append((name, "as written", settings, lines))
        if not args.no_variants:
            for action, what, new in variants(lines):
                # the same change of the same keyword of the same action is tried once only
                key = (action, what.split(": ", 1)[1])
                if args.all_variants or key not in seen:
                    seen.add(key)
                    runs.append((name, what, settings, new))
    if args.run:
        os.makedirs(args.out, exist_ok=True)
        args.trajectory = os.path.abspath(os.path.join(args.out, "traj.xyz"))
        make_trajectory(args.out, 3000)
    print(f"{len(examples)} examples ({skipped} skipped), {len(runs)} runs", flush=True)

    # the files the examples read are taken from a copy of regtest/ and new-manual/extras, which the
    # examples can modify without touching the source tree
    args.data = tempfile.mkdtemp(prefix="auditdata")
    shutil.copytree(os.path.join(ROOT, "regtest"), os.path.join(args.data, "regtest"), symlinks=True,
                    ignore=shutil.ignore_patterns("tmp"))
    shutil.copytree(os.path.join(ROOT, "new-manual", "extras"), os.path.join(args.data, "extras"))

    found = []
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = {pool.submit(run_one, r[3], r[2], args): r for r in runs}
        for n, fut in enumerate(concurrent.futures.as_completed(futures)):
            if n % 500 == 499:
                print(f"{n + 1} runs done, {len(found)} problems so far", flush=True)
            name, what, settings, lines = futures[fut]
            kind, output = fut.result()
            if kind:
                found.append((kind, name, what, lines, output))

    shutil.rmtree(args.data, ignore_errors=True)
    os.makedirs(args.out, exist_ok=True)
    found.sort(key=lambda x: (x[0], x[1], x[2]))
    for n, (kind, name, what, lines, output) in enumerate(found):
        with open(os.path.join(args.out, f"{n:04d}-{kind}.txt"), "w") as f:
            f.write(f"# {kind} {name} ({what})\n" + "\n".join(lines) + "\n# ---- output\n" + output)
        message = ""
        for line in output.splitlines():
            if any(p.search(line) for _, p in SIGNATURES):
                message = line.strip()[:200]
                break
        print(f"{kind:9s} {name} ({what}): {message}")
    print(f"{len(found)} problems in {len(runs)} runs; details in {args.out}/")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
