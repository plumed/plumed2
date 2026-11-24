#!/usr/bin/env python3
"""Audit installed PLUMED ELF imports against its declared Linux ABI baseline."""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", type=Path, required=True)
args = parser.parse_args()
if sys.platform != "linux":
    print("ELF ABI audit applies only on Linux")
    raise SystemExit(0)
prefix = args.prefix.resolve()
records = [json.loads(p.read_text()) for p in (prefix / "conda-meta").glob("*.json")]
packages = {r["name"] for r in records}
count = 0
libraries = set((prefix / "lib").glob("libplumed*.so*"))
libraries.update((prefix / "lib").glob("*Coordination.so"))
libraries.update((prefix / "lib").glob("*OpenACC.so"))
for path in sorted(libraries):
    if path.is_symlink() or not path.is_file():
        continue
    versions = subprocess.check_output(
        ["readelf", "--version-info", str(path)], text=True
    )
    imports = {
        tuple(map(int, v.split(".")))
        for v in re.findall(r"\bGLIBC_([0-9.]+)", versions)
    }
    maximum = max(imports, default=(0,))
    assert maximum <= (2, 28), (
        f"{path.name} imports GLIBC {maximum}, exceeds declared 2.28"
    )
    dynamic = subprocess.check_output(["readelf", "-d", str(path)], text=True)
    needed = re.findall(r"\(NEEDED\).*?\[(.*?)\]", dynamic)
    if any("openblas" in lib for lib in needed):
        assert "libopenblas" in packages, (
            f"{path.name}: missing actual OpenBLAS provider"
        )
    assert not any(lib in ("libblas.so.3", "liblapack.so.3") for lib in needed), needed
    print(
        f"{path.name}: maximum GLIBC {'.'.join(map(str, maximum))}; providers checked"
    )
    count += 1
assert count, "No PLUMED ELF libraries found"
