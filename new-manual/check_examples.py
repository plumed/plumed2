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


# MPI (OpenMPI) sometimes fails to start when many examples are run at the same time: such a failure
# says nothing about the example and does not count
MPI_STARTUP = ("orte_init failed", "ORTE_ERROR_LOG", "PMIX_ERROR", "MPI_Init")


def first_error(filename):
    with open(filename, errors="replace") as f:
        lines = f.read().splitlines()
    for i, line in enumerate(lines):
        if any(m in line for m in MPI_STARTUP):
            return "MPI", " ".join(lines[i : i + 3])[:400]
        if "ERROR" in line or "assertion failed" in line or "Segmentation" in line:
            return "PLUMED", " ".join(lines[i : i + 3])[:400]
    return None, ""


def stderr_files(page):
    # PlumedToHTML writes the stderr of each example as <page>.md_working_<n>.dat.plumed.stderr.md;
    # the pages of modules and plugins are module_<name>.md and plugins_<name>.md
    name = page.removesuffix(".md")
    files = glob.glob(f"docs/{name}.md_working_*.plumed.stderr.md")
    files += glob.glob(f"docs/*_{name}.md_working_*.plumed.stderr.md")
    return sorted(files)


broken = read("broken-examples.txt")
allowed = read("broken-examples.allowed")

failed = False
for page, n in sorted(broken.items()):
    if n <= allowed.get(page, 0):
        continue
    errors = [(f,) + first_error(f) for f in stderr_files(page)]
    mpi = sum(1 for e in errors if e[1] == "MPI")
    real = n - mpi > allowed.get(page, 0)
    failed = failed or real
    print(
        f"{'ERROR' if real else 'WARNING'}: page {page} has {n} failing examples,"
        f" {allowed.get(page, 0)} allowed"
        + (f" ({mpi} because MPI did not start, not counted)" if mpi else "")
    )
    for f, kind, msg in errors:
        if kind:
            print(f"  {f}: {msg}")
for page, n in sorted(allowed.items()):
    if broken.get(page, 0) < n:
        print(
            f"NOTE: page {page} has {broken.get(page, 0)} failing examples, {n} allowed:"
            " please update broken-examples.allowed"
        )

print(f"{sum(broken.values())} failing examples in {len(broken)} pages")
if failed:
    print(
        "Fix the examples above (or, for an example that cannot run here, use the #SETTINGS"
        " options of PlumedToHTML), and do not add them to broken-examples.allowed"
    )
    sys.exit(1)
