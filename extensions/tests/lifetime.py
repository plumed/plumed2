#!/usr/bin/env python3
"""Exercise native extension teardown through the runtime wrapper."""

import argparse
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path


def check_lifetime(base):
    source = Path(__file__).resolve().parents[2]
    base = base.resolve()
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    kernel = base / "lib" / ("libplumedKernel" + suffix)
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    common = ["-std=c++17", "-I" + str(base / "include")]
    module_flags = ["-dynamiclib"] if sys.platform == "darwin" else ["-shared", "-fPIC"]
    # GNU-unique symbols would independently mark the fixture NODELETE on ELF.
    probe = subprocess.run(
        [*compiler, "-fno-gnu-unique", "-x", "c++", "-c", "-o", os.devnull, "-"],
        input="",
        text=True,
        capture_output=True,
        check=False,
    )
    if probe.returncode == 0:
        module_flags.append("-fno-gnu-unique")
    with tempfile.TemporaryDirectory(prefix="plumed-module-lifetime-") as scratch:
        root = Path(scratch)
        libraries = {}
        for name in ("lifetime_dep", "lifetime_root", "zbroken"):
            library = root / (name + suffix)
            libraries[name] = library
            subprocess.run(
                [
                    sys.executable,
                    str(source / "extensions/register.py"),
                    "--plumed",
                    str(base / "bin/plumed"),
                    "--module",
                    name,
                    "--library",
                    str(library),
                    "--output",
                    str(root / name),
                    "--requires",
                    *(["lifetime_dep"] if name == "lifetime_root" else []),
                ],
                check=True,
                timeout=30,
            )
            subprocess.run(
                [
                    *compiler,
                    *common,
                    *module_flags,
                    "-I" + str(base / "include/plumed"),
                    '-DLIFETIME_NAME="' + name + '"',
                    *(["-DLIFETIME_ROOT"] if name == "lifetime_root" else []),
                    str(root / name / "module_identity.cpp"),
                    str(source / "extensions/tests/lifetime_module.cpp"),
                    *(
                        [str(libraries["lifetime_dep"])]
                        if name == "lifetime_root"
                        else []
                    ),
                    str(kernel),
                    "-Wl,-rpath," + str(base / "lib"),
                    "-o",
                    str(library),
                ],
                check=True,
                timeout=120,
            )
        program = root / "lifetime"
        subprocess.run(
            [
                *compiler,
                *common,
                str(source / "extensions/tests/lifetime.cpp"),
                "-pthread",
                *([] if sys.platform == "darwin" else ["-ldl"]),
                "-o",
                str(program),
            ],
            check=True,
            timeout=120,
        )
        catalog = root / "catalog"
        catalog.mkdir()
        for name in ("lifetime_dep", "lifetime_root"):
            (catalog / (name + ".plumed")).write_text(
                (root / name / (name + ".plumed")).read_text()
            )
        direct = root / "lifetime-cli"
        subprocess.run(
            [
                *compiler,
                *common,
                "-DLIFETIME_DIRECT_CLI",
                str(source / "extensions/tests/lifetime.cpp"),
                str(kernel),
                "-Wl,-rpath," + str(base / "lib"),
                "-pthread",
                *([] if sys.platform == "darwin" else ["-ldl"]),
                "-o",
                str(direct),
            ],
            check=True,
            timeout=120,
        )
        for mode, cycles in (
            ("overlap", 3),
            ("threads", 1),
            ("disabled", 1),
            ("cli", 1),
            ("failure", 2),
        ):
            events = root / (mode + ".events")
            extra = []
            if mode == "failure":
                (catalog / "zbroken.plumed").write_text(
                    (root / "zbroken/zbroken.plumed").read_text()
                )
                hidden = root / "zbroken.hidden"
                libraries["zbroken"].rename(hidden)
                extra = [str(hidden), str(libraries["zbroken"])]
            subprocess.run(
                [
                    str(direct if mode == "cli" else program),
                    mode,
                    str(kernel),
                    str(libraries["lifetime_root"]),
                    *extra,
                ],
                env=dict(
                    os.environ,
                    PLUMED_AUTOLOAD="1",
                    PLUMED_PLUGIN_PATH=str(catalog),
                    PLUMED_TEST_LIFETIME_EVENTS=str(events),
                ),
                check=True,
                timeout=120,
            )
            # Each dependent must finish destruction before its dependency.
            expected = ["lifetime_root", "lifetime_dep"] * cycles
            if mode == "failure":
                expected.insert(2, "zbroken")
            assert events.read_text().splitlines() == expected, events.read_text()
    print("Autoload lifetime, rollback and native dependency order passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    check_lifetime(parser.parse_args().base)
