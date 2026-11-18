#!/usr/bin/env python3
"""Check kernel/extension build-feature vectors without cross-module leakage."""

import argparse
import importlib
import json
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--base", type=Path, required=True)
args = parser.parse_args()
base = args.base.resolve()
program = str(base / "bin/plumed")
source = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(source / "extensions"))
module_symbol = importlib.import_module("symbols").module_symbol

env = dict(os.environ, PLUMED_AUTOLOAD="0")
soext = subprocess.check_output(
    [program, "info", "--soext"], text=True, env=env
).strip()
compiler = shlex.split(os.environ.get("CXX", "c++"))

with tempfile.TemporaryDirectory(prefix="plumed-module-features-") as scratch:
    root = Path(scratch)
    env["PLUMED_PLUGIN_PATH"] = str(root)
    for name, feature, dependencies in (
        ("feature-a", "TEST_ALPHA", []),
        ("feature_a", "TEST_UNDERSCORE", []),
        ("feature_b", "TEST_BETA", ["feature-a"]),
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
                    f'extern "C" const char* {module_symbol(dependencies[0])}();\n'
                    f'extern "C" const void* feature_dependency() {{ return reinterpret_cast<const void*>(&{module_symbol(dependencies[0])}); }}\n'
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
                *(
                    [str(source / "extensions/tests/metadata.cpp")]
                    if name == "feature-a"
                    else []
                ),
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
        query("feature-a", feature, 1)
    # A hyphenated name and its underscore spelling have distinct native symbols.
    query("feature_a", "test_underscore", 1)
    isolated = subprocess.run(
        [
            program,
            "--load-module=feature-a",
            "--load-module=feature_a",
            "modules",
            "--module",
            "feature_a",
            "--has-feature",
            "test_underscore",
            "-q",
        ],
        env=env,
        check=False,
        timeout=30,
    )
    assert isolated.returncode == 0
    query("feature-a", "test_underscore", 1)
    query("feature-a", "test_alpha", 0)
    query("feature_b", "test_beta", 0)
    query("feature_b", "test_alpha", 1)
    query("feature-a", "test_beta", 1)
    query("tools", "test_alpha", 1)
    query("missing", "test_alpha", 1)
    query("feature-a", "test_alpha", 1, load=False)
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
    # The same compiled modules must also be discoverable through eager loading,
    # including their native dependency and distinct hyphen/underscore names.
    subprocess.run(
        [
            program,
            "modules",
            "--has",
            "feature-a",
            "feature_a",
            "feature_b",
            "-q",
        ],
        env=dict(env, PLUMED_AUTOLOAD="1"),
        check=True,
        timeout=30,
    )

    # All inquiry clients must reflect the selected runtime context, including
    # native dependencies, with logical names rather than image-address keys.
    def metadata(*arguments, load=True, eager=False):
        return subprocess.run(
            [
                program,
                "--no-mpi",
                *(["--load-module=feature_b"] if load else []),
                *arguments,
            ],
            env=dict(env, PLUMED_AUTOLOAD="1" if eager else "0"),
            check=True,
            text=True,
            capture_output=True,
            timeout=30,
        )

    before = metadata("gentemplate", "--list", load=False).stderr.split()
    assert "METADATA_ACTION" not in before
    for eager in (False, True):
        listed = metadata(
            "gentemplate", "--list", load=not eager, eager=eager
        ).stderr.split()
        assert "METADATA_ACTION" in listed, listed
        assert not any(":" in name for name in listed), listed
        manual = metadata(
            "manual",
            "--action",
            "METADATA_ACTION",
            "--vim",
            load=not eager,
            eager=eager,
        )
        assert manual.stdout.startswith("METADATA_ACTION,")
        assert "metadata_tool" in manual.stderr.split()
        assert not any(name.startswith("0x") for name in manual.stderr.split())
        assert (
            "[VALUE]"
            in metadata(
                "gentemplate",
                "--action",
                "METADATA_ACTION",
                "--include-optional",
                load=not eager,
                eager=eager,
            ).stdout
        )
        document = json.loads(metadata("gen_json", load=not eager, eager=eager).stdout)
        assert "METADATA_ACTION" in document
        assert "metadata_tool" in document["cltools"]
        metadata("metadata_tool", "--probe", load=not eager, eager=eager)
        destination = root / ("vim-eager" if eager else "vim-explicit")
        metadata("gen_vim", "--outdir", str(destination), load=not eager, eager=eager)
        assert (
            'b:plumedDictionary["METADATA_ACTION"]'
            in (destination / "syntax/plumed.vim").read_text()
        )
        assert (destination / "help/METADATA_ACTION.txt").is_file()
    # Real Bash completion must carry explicit LOAD options into its metadata
    # query, and expose the modules command and its presence/feature switches.
    completion_env = dict(env, PATH=str(base / "bin") + os.pathsep + os.environ["PATH"])
    for words, expected in (
        (["plumed", "--"], "--load-module"),
        (["plumed", "modules", "--"], "--has-feature"),
        (["plumed", "--load-module", "feature_b", "metadata_tool", "--"], "--probe"),
        (["plumed", "--load-module=feature_b", "metadata_tool", "--"], "--probe"),
        (["plumed", "driver", "--load", str(root / "missing.so"), "--"], "--plumed"),
        (["plumed", "driver", "--load=" + str(root / "missing.so"), "--"], "--plumed"),
        (["plumed", "driver", "--load-module", "missing-module", "--"], "--plumed"),
        (
            [
                "plumed",
                "--load-module=feature_b",
                "driver",
                "--load-module=missing-module",
                "--",
            ],
            "--plumed",
        ),
    ):
        script = 'source "$1"\nshift\nCOMP_WORDS=("$@")\nCOMP_CWORD=$((${#COMP_WORDS[@]}-1))\n_plumed\nprintf "%s\\n" "${COMPREPLY[@]}"'
        result = subprocess.run(
            [
                "bash",
                "-c",
                script,
                "completion-test",
                str(source / "src/lib/bash.completion"),
                *words,
            ],
            env=completion_env,
            check=True,
            text=True,
            capture_output=True,
            timeout=30,
        )
        assert expected in result.stdout.split(), (words, result.stdout, result.stderr)
print("Per-module kernel/extension feature isolation passed")
