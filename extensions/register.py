#!/usr/bin/env python3
"""Generate registration artifacts for an existing plugin's own build system."""

import argparse
import os
import re
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--plumed", required=True)
parser.add_argument("--module", required=True)
parser.add_argument(
    "--library",
    required=True,
    help="library path relative to the installed registration directory",
)
parser.add_argument("--requires", nargs="*", default=[])
parser.add_argument("--variant", default="default")
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()
for name in [args.module, *args.requires]:
    if not re.fullmatch(r"[a-z][a-z0-9_]*", name):
        parser.error("module names must be lowercase C identifiers")
for text in [args.library, args.variant]:
    if "\n" in text or "\r" in text:
        parser.error("registration values must be single lines")
build_id = subprocess.check_output(
    [args.plumed, "--module-build-id"],
    text=True,
    env=dict(os.environ, PLUMED_AUTOLOAD="0"),
).strip()
args.output.mkdir(parents=True, exist_ok=True)
(args.output / "module_identity.cpp").write_text(
    '#include "config/Config.h"\n'
    + f'extern "C" const char* plumed_module_{args.module}() {{ return "{build_id}"; }}\n'
    + f'extern "C" const void* plumed_module_{args.module}_kernel() {{ return reinterpret_cast<const void*>(&PLMD::config::getModuleBuildID); }}\n'
)
(args.output / f"{args.module}.plumed").write_text(
    f"module {args.module}\nbuild {build_id}\nvariant {args.variant}\nlibrary {args.library}\n"
    + "requires "
    + " ".join(args.requires)
    + "\n"
)
