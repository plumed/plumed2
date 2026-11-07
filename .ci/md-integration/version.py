#!/usr/bin/env python3
"""Read an engine's pinned source version, locally or as GitHub step outputs."""

import argparse
import json
import os
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("engine", choices=("gromacs", "lammps", "amber"))
    parser.add_argument("--github-output", action="store_true")
    args = parser.parse_args()
    versions = json.loads(Path(__file__).with_name("versions.json").read_text())
    entry = versions[args.engine]
    digest = "sha256" if args.engine == "amber" else "ref"
    length = 64 if digest == "sha256" else 40
    if not re.fullmatch(r"[0-9a-f]{%d}" % length, entry[digest]):
        parser.error(f"{args.engine}: {digest} must pin an exact source digest")
    if any(
        not isinstance(value, str) or "\n" in value or "\r" in value
        for value in entry.values()
    ):
        parser.error("Version fields must be single-line strings")
    print(json.dumps(entry, indent=2))
    if args.github_output:
        with open(os.environ["GITHUB_OUTPUT"], "a") as stream:
            for key, value in entry.items():
                stream.write(f"{key}={value}\n")


if __name__ == "__main__":
    main()
