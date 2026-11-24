"""Ensure all distributed GROMACS patch originals match the pinned source."""

import subprocess
import sys
from pathlib import Path

source = Path(sys.argv[1]).resolve()
patch = Path(__file__).resolve().parents[3] / "patches/gromacs-2026.4.diff"
originals = sorted(patch.rglob("*.preplumed"))
if not originals:
    raise RuntimeError("No distributed GROMACS patch originals found")
for original in originals:
    relative = original.relative_to(patch).as_posix().removesuffix(".preplumed")
    expected = subprocess.run(
        ["git", "-C", str(source), "show", f"HEAD:{relative}"],
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    if original.read_bytes() != expected:
        raise RuntimeError(f"Patch original differs from pinned GROMACS: {relative}")
print(f"All {len(originals)} distributed originals match pinned GROMACS")
