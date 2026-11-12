#!/usr/bin/env python3
"""Check kernel/extension build-feature vectors without cross-module leakage."""

import argparse
import os
import shlex
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--base", type=Path, required=True)
args = parser.parse_args()
base = args.base.resolve()
program = str(base / "bin/plumed")
source = Path(__file__).resolve().parents[2]
env = dict(os.environ, PLUMED_AUTOLOAD="0")
soext = subprocess.check_output(
    [program, "info", "--soext"], text=True, env=env
).strip()
compiler = shlex.split(os.environ.get("CXX", "c++"))

with tempfile.TemporaryDirectory(prefix="plumed-module-features-") as scratch:
    root = Path(scratch)
    env["PLUMED_PLUGIN_PATH"] = str(root)
    for name, feature, dependencies in (
        ("feature_a", "TEST_ALPHA", []),
        ("feature_b", "TEST_BETA", ["feature_a"]),
    ):
        subprocess.run(
            [
                os.sys.executable,
                str(source / "extensions/register.py"),
                "--plumed",
                program,
                "--module",
                name,
                "--library",
                str(root / (name + "." + soext)),
                "--features",
                feature,
                "--requires",
                *dependencies,
                "--output",
                str(root / name),
            ],
            check=True,
            env=env,
            timeout=30,
        )
        identity = root / name / "module_identity.cpp"
        if dependencies:
            with identity.open("a") as handle:
                handle.write(
                    'extern "C" const char* plumed_module_feature_a();\nextern "C" const void* feature_dependency() { return reinterpret_cast<const void*>(&plumed_module_feature_a); }\n'
                )
        flags = ["-dynamiclib"] if soext == "dylib" else ["-shared", "-fPIC"]
        subprocess.run(
            [
                *compiler,
                *flags,
                "-std=c++17",
                "-D__PLUMED_HAS_" + feature + "=1",
                "-I" + str(base / "include/plumed"),
                str(root / name / "module_identity.cpp"),
                *[
                    str(root / (dependency + "." + soext))
                    for dependency in dependencies
                ],
                str(base / "lib" / ("libplumedKernel." + soext)),
                "-Wl,-rpath," + str(base / "lib"),
                "-o",
                str(root / (name + "." + soext)),
            ],
            check=True,
            timeout=120,
        )
        (root / (name + ".plumed")).write_text(
            (root / name / (name + ".plumed")).read_text()
        )

    def query(module, feature, expected, load=True):
        command = [program] + (["--load-module=feature_b"] if load else [])
        result = subprocess.run(
            [*command, "modules", "--module", module, "--has-feature", feature, "-q"],
            env=env,
            check=False,
            timeout=30,
        )
        assert result.returncode == expected, (module, feature, result.returncode)

    # Builtin queries return exactly the kernel vector.
    baseline = subprocess.check_output(
        [program, "modules", "--module", "tools", "--list"], env=env, text=True
    )
    aggregate = subprocess.check_output(
        [program, "modules", "--list"], env=env, text=True
    )
    assert "tools" in aggregate.split()
    for feature in baseline.split():
        query("tools", feature, 0)
        query("feature_a", feature, 1)
    query("feature_a", "test_alpha", 0)
    query("feature_b", "test_beta", 0)
    query("feature_b", "test_alpha", 1)
    query("feature_a", "test_beta", 1)
    query("tools", "test_alpha", 1)
    query("missing", "test_alpha", 1)
    query("feature_a", "test_alpha", 1, load=False)
    subprocess.run(
        [
            program,
            "--load-module=feature_b",
            "modules",
            "--has-feature",
            "test_alpha",
            "test_beta",
            "-q",
        ],
        env=env,
        check=True,
        timeout=30,
    )
print("Per-module kernel/extension feature isolation passed")
