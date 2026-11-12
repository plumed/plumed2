#!/usr/bin/env python3
"""Check CMake selection of hyphenated modules after upstream integration."""

import os
import subprocess
import tempfile
from pathlib import Path

source = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="plumed-cmake-modules-") as directory:
    for expression, present, absent in (
        ("reset:+isdb-pytorch", {"isdb", "isdb-pytorch"}, set()),
        ("all:-isdb-pytorch", {"isdb", "ves"}, {"isdb-pytorch"}),
        ("reset+ves-isdb", {"ves"}, {"isdb", "isdb-pytorch"}),
    ):
        result = subprocess.run(
            [
                os.environ.get("CMAKE", "cmake"),
                "-S",
                str(source),
                "-B",
                directory,
                "-DPLUMED_MODULES=" + expression,
            ],
            capture_output=True,
            text=True,
            check=False,
            timeout=120,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        selected = next(
            line.removeprefix("-- PLUMED modules: ").split()
            for line in result.stdout.splitlines()
            if line.startswith("-- PLUMED modules: ")
        )
        assert present <= set(selected) and not absent & set(selected), (
            expression,
            selected,
        )
print("CMake hyphenated and legacy module selection passed")
