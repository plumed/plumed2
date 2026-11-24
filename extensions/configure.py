#!/usr/bin/env python3
"""Configure extension-only builds against an installed PLUMED base."""

import argparse
import os
import shlex
import shutil
import subprocess
from pathlib import Path

from inventory import generate_catalog

root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser()
parser.add_argument("--with-extension-base", required=True, type=Path)
parser.add_argument("--enable-extension-modules", required=True)
parser.add_argument("--prefix", type=Path)
# Modules with "toolchain: cuda" (as for the native configure): nvcc by default,
# hipcc with --enable-hip. The CUDA runtime is shared by default, as for the
# conda packages (one package per CUDA version on the same base).
parser.add_argument("--enable-hip", action="store_true")
parser.add_argument("--with-cudart", choices=("shared", "static"), default="shared")
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
catalog = generate_catalog()
modules = args.enable_extension_modules.replace(":", " ").split()
if modules == ["all"]:
    # Modules compiled with another toolchain are built only when named.
    modules = [
        name
        for name, module in catalog["modules"].items()
        if module["type"] == "default-off"
        and module["toolchain"] == "cxx"
        and name not in builtin
    ]
if not modules:
    parser.error("select at least one extension module")
for name in modules:
    if (
        not name.replace("_", "").replace("-", "").isalnum()
        or name not in catalog["modules"]
    ):
        parser.error(f"unknown module: {name}")
    if name in builtin:
        parser.error(f"{name} is already built into the selected base")

def closure(name, seen):
    if name in seen or name in builtin or name not in catalog["modules"]:
        return seen
    seen.add(name)
    for dependency in catalog["modules"][name]["source_use"]:
        closure(dependency, seen)
    return seen


selection = set()
for name in modules:
    closure(name, selection)
cuda = []
if any(catalog["modules"][name]["toolchain"] == "cuda" for name in selection):
    backend = "hipcc" if args.enable_hip else "nvcc"
    variable = "HIPCC" if args.enable_hip else "NVCC"
    command = shlex.split(
        os.environ.get(variable)
        or ("" if args.enable_hip else os.environ.get("CUDACXX", ""))
        or backend
    )
    if not shutil.which(command[0]):
        parser.error(
            f"{variable}={command[0]} not found: it is needed by the modules with toolchain cuda"
        )
    cuda = [
        f"EXTENSION_CUDA_BACKEND={backend}",
        f"EXTENSION_CUDA_COMPILER={shlex.join(command)}",
        f"EXTENSION_CUDART={args.with_cudart}",
        f"EXTENSION_NVCCFLAGS={os.environ.get('NVCCFLAGS', '')}",
        f"EXTENSION_HIP_ARCHITECTURES={os.environ.get('HIP_ARCHITECTURES', '')}",
        f"EXTENSION_CUDA_LIBDIR={os.environ.get('CUDA_LIBDIR', '')}",
    ]
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
            *cuda,
        ]
    )
    + "\n"
)
print(f"Configured extensions {', '.join(modules)} against base {build_id}")
print("Run make, then make install. No kernel will be rebuilt or installed.")
