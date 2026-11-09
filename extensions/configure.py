#!/usr/bin/env python3
"""Configure extension-only builds against an installed PLUMED base."""

import argparse
import json
import os
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser()
parser.add_argument("--with-extension-base", required=True, type=Path)
parser.add_argument("--enable-extension-modules", required=True)
parser.add_argument("--prefix", type=Path)
args = parser.parse_args()
base = args.with_extension_base.resolve()
prefix = (args.prefix or base).resolve()
environment = dict(os.environ, PLUMED_AUTOLOAD="0")
plumed = base / "bin/plumed"
config = subprocess.check_output(
    [str(plumed), "config", "makefile_conf"], env=environment, text=True
)
builtin = set(
    subprocess.check_output(
        [str(plumed), "--list-modules"], env=environment, text=True
    ).split()
)
build_id = subprocess.check_output(
    [str(plumed), "--module-build-id"], env=environment, text=True
).strip()
catalog = json.loads((root / "extensions/catalog.json").read_text())
features = set(re.findall(r"-D__PLUMED_HAS_([A-Z0-9_]+)(?:=1)?", config))
modules = args.enable_extension_modules.replace(":", " ").split()
if modules == ["all"]:
    modules = [
        p.parent.name
        for p in sorted((root / "src").glob("*/module.type"))
        if p.read_text().strip() == "default-off"
        and set(catalog["modules"][p.parent.name]["required_features"]) <= features
    ]
if not modules:
    parser.error("select at least one extension module")
for name in modules:
    if (
        not name.replace("_", "").isalnum()
        or not (root / "src" / name / "module.type").is_file()
    ):
        parser.error(f"unknown module: {name}")
    missing = set(catalog["modules"][name]["required_features"]) - features
    if missing:
        parser.error(
            f"{name} needs base feature macros: {sorted(missing)}; dependency-enabled extension builds are not yet configured"
        )
    if name in builtin:
        parser.error(f"{name} is already built into the selected base")
# GNU Make paths are intentionally constrained; shell commands quote all paths.
for path in (base, prefix, root):
    if any(c in str(path) for c in " \t\n#$"):
        parser.error("extension build paths must not contain whitespace, #, or $")
(root / "Makefile.conf").write_text(
    config
    + "\n"
    + "\n".join(
        [
            f"EXTENSION_BASE={base}",
            f"EXTENSION_PREFIX={prefix}",
            f"EXTENSION_MODULES={' '.join(modules)}",
            f"EXTENSION_BUILD_ID={build_id}",
            f"EXTENSION_PYTHON={os.sys.executable}",
        ]
    )
    + "\n"
)
print(f"Configured extensions {', '.join(modules)} against base {build_id}")
print("Run make, then make install. No kernel will be rebuilt or installed.")
