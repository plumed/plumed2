#!/usr/bin/env python3
"""Check real configure selection for hyphenated names and legacy expressions."""

import shutil
import subprocess
import tempfile
from pathlib import Path

source = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="plumed-configure-modules-") as scratch:
    root = Path(scratch)
    shutil.copy2(source / "configure", root / "configure")
    for declaration in (source / "src").glob("*/module.type"):
        target = root / declaration.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(declaration, target)
    for expression, enabled, disabled in (
        ("reset:+isdb-pytorch", {"isdb-pytorch"}, set()),
        ("all:-isdb-pytorch", set(), {"isdb-pytorch"}),
        ("reset+ves-isdb", {"ves"}, {"isdb"}),
        ("none:isdb-pytorch:+pytorch", {"isdb-pytorch", "pytorch"}, {"isdb"}),
    ):
        # Fail deliberately after selection, before any compiler/dependency probes.
        result = subprocess.run(
            ["bash", "configure", "--enable-modules=" + expression + ":missing_probe"],
            cwd=root,
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
        assert result.returncode != 0
        assert "module missing_probe which does not exist" in result.stderr, (
            result.stderr
        )
        for name in enabled:
            assert (root / "src" / (name + ".on")).exists(), (expression, name)
            assert not (root / "src" / (name + ".off")).exists(), (expression, name)
        for name in disabled:
            assert (root / "src" / (name + ".off")).exists(), (expression, name)
            assert not (root / "src" / (name + ".on")).exists(), (expression, name)
print("Hyphenated and legacy native configure module selection passed")
