#!/usr/bin/env python3
"""Record public thread-MPI header changes using a local GROMACS git clone.

This is an audit aid, not a declaration that every discovered API is safe to
publish. Review changed signatures and implementations before updating guards.
"""

import argparse
import difflib
import json
from pathlib import Path
import re
import subprocess

HEADERS = (
    "src/external/thread_mpi/include/thread_mpi/tmpi.h",
    "src/external/thread_mpi/include/thread_mpi/mpi_bindings.h",
)


def git(checkout, *arguments):
    return subprocess.check_output(["git", "-C", str(checkout), *arguments], text=True)


def release_key(tag):
    parts = tag.removeprefix("v").split(".")
    return int(parts[0]), int(parts[1]) if len(parts) > 1 else 0


def public_lines(text):
    # Comments, whitespace and const formatting alone are not declarations.
    # Preserve const tokens and all non-comment source, including typedefs and
    # macros: function names alone would miss the Error_string ABI change.
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//[^\n]*", "", text)
    return [line for line in (" ".join(x.split()) for x in text.splitlines()) if line]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gromacs", type=Path)
    parser.add_argument(
        "--through", default="v2026.4", help="last stable release to inspect"
    )
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("gromacs-history.json")
    )
    parser.add_argument("--check", action="store_true", help="compare without writing")
    args = parser.parse_args()
    tags = git(args.gromacs, "tag", "--list", "v20*").splitlines()
    tags = sorted(
        (
            tag
            for tag in tags
            if re.fullmatch(r"v20\d\d(?:\.\d+)?", tag)
            and (2020, 0) <= release_key(tag) <= release_key(args.through)
        ),
        key=release_key,
    )
    if not tags or tags[0] != "v2020" or args.through not in tags:
        parser.error(
            "local clone must contain v2020 and the requested final stable release"
        )
    previous = None
    changes = []
    for tag in tags:
        source = {
            path: public_lines(git(args.gromacs, "show", f"{tag}:{path}"))
            for path in HEADERS
        }
        if source != previous:
            row = {
                "first_observed_release": tag,
                "commit": git(args.gromacs, "rev-parse", f"{tag}^{{commit}}").strip(),
            }
            bindings = git(args.gromacs, "show", f"{tag}:{HEADERS[1]}")
            aliases = dict(
                re.findall(r"^#define\s+(MPI_\w+)\s+(tMPI_\w+)\b", bindings, re.M)
            )
            row["mpi_function_aliases"] = aliases
            row["header_changes"] = {
                path: list(
                    difflib.unified_diff(
                        previous[path] if previous else [],
                        source[path],
                        fromfile="previous inspected stable release",
                        tofile=tag,
                        lineterm="",
                    )
                )
                for path in HEADERS
                if previous is None or source[path] != previous[path]
            }
            changes.append(row)
        previous = source
    report = {
        "schema_version": 1,
        "baseline": "v2020",
        "through": args.through,
        "note": "Observed declarations/bindings only; baseline APIs may predate 2020. Review semantics separately.",
        "inspected_releases": tags,
        "changes": changes,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit(
                "thread-MPI history differs; regenerate and review API compatibility guards"
            )
        print(f"thread-MPI history unchanged across {len(tags)} stable releases")
    else:
        args.output.write_text(rendered)
        print(
            f"recorded {len(changes)} public-header revisions across {len(tags)} stable releases"
        )


if __name__ == "__main__":
    main()
