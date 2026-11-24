#!/usr/bin/env python3
"""Verify EMMIVOX belongs only to the separately loaded ISDB PyTorch module."""

import argparse
import json
import os
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", type=Path, required=True)
args = parser.parse_args()
prefix = args.prefix.resolve()
program = str(prefix / "bin/plumed")
report = json.loads((prefix / "lib/plumed/plugins/isdb-pytorch.json").read_text())
assert report["optional_features"]["LIBTORCH"], report


def run(*arguments, load=False, autoload=False, expected=0):
    command = [program]
    if load:
        command.append("--load-module=isdb-pytorch")
    result = subprocess.run(
        [*command, *arguments],
        env=dict(os.environ, PLUMED_AUTOLOAD="1" if autoload else "0"),
        text=True,
        capture_output=True,
        check=False,
        timeout=60,
    )
    assert result.returncode == expected, (
        command,
        arguments,
        result.stdout,
        result.stderr,
    )
    return result.stdout


# The kernel keeps ordinary ISDB functionality but no EMMIVOX/Torch dependency.
run("modules", "--has", "isdb", "-q")
run("modules", "--has", "isdb-pytorch", "-q", expected=1)
run("manual", "--action", "EMMI")
run("manual", "--action", "EMMIVOX", expected=1)
for mode in ({"load": True}, {"autoload": True}):
    run("modules", "--has", "isdb-pytorch", "-q", **mode)
    run("manual", "--action", "EMMIVOX", **mode)
    run(
        "modules", "--module", "isdb-pytorch", "--has-feature", "libtorch", "-q", **mode
    )
    run(
        "modules",
        "--module",
        "isdb",
        "--has-feature",
        "libtorch",
        "-q",
        expected=1,
        **mode,
    )
    run("config", "has", "libtorch", "-q", expected=1, **mode)
# EMMIVOX does not introduce a hard dependency on the PLUMED pytorch module.
run("modules", "--has", "pytorch", "-q", load=True, expected=1)
print(
    "ISDB PyTorch action ownership, automatic/explicit loading and feature isolation passed"
)
