#!/usr/bin/env python3
"""Fail when an example of the manual that used to work does not work any more.

Run after "make doc" (build_manual.py writes broken-examples.txt). Every page with failing examples
must be listed in broken-examples.allowed with at most that number of failures. When an example is
fixed, lower or remove its line in broken-examples.allowed, otherwise this script asks you to.
"""

import glob
import sys


def read(filename):
    pages = {}
    with open(filename) as f:
        for line in f:
            line = line.split("#")[0].split()
            if line:
                pages[line[0]] = int(line[1])
    return pages


def first_error(filename):
    with open(filename, errors="replace") as f:
        lines = f.read().splitlines()
    for i, line in enumerate(lines):
        if "ERROR" in line or "assertion failed" in line or "Segmentation" in line:
            return " ".join(lines[i : i + 3])[:400]
    return ""


broken = read("broken-examples.txt")
allowed = read("broken-examples.allowed")

new = {p: n for p, n in broken.items() if n > allowed.get(p, 0)}
fixed = {p: n for p, n in allowed.items() if broken.get(p, 0) < n}

for page, n in sorted(new.items()):
    print(f"ERROR: page {page} has {n} failing examples, {allowed.get(page, 0)} allowed")
    for err in sorted(glob.glob(f"docs/*{page.removesuffix('.md')}.md_working_*.plumed.stderr")):
        msg = first_error(err)
        if msg:
            print(f"  {err}: {msg}")
for page, n in sorted(fixed.items()):
    print(f"NOTE: page {page} has {broken.get(page, 0)} failing examples, {n} allowed:"
          " please update broken-examples.allowed")

print(f"{sum(broken.values())} failing examples in {len(broken)} pages")
if new:
    print("Fix the examples above (or, for an example that cannot run here, use the #SETTINGS"
          " options of PlumedToHTML), and do not add them to broken-examples.allowed")
    sys.exit(1)
