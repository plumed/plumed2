#!/usr/bin/env python3
"""Exercise module packaging declarations, generated policy and build backends."""

import argparse
import copy
import hashlib
import importlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import yaml

source = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(source / "extensions"))
config_module = importlib.import_module("conda_config")
build_environment = config_module.build_environment
merge_policies = config_module.merge_policies
regression_selection = config_module.regression_selection
run_build = config_module.run_build
skip_selector = config_module.skip_selector
supported = config_module.supported
write_variants = config_module.write_variants
policy_module = importlib.import_module("conda_policy")
read_policy = policy_module.read_policy
validate_backend = policy_module.validate_backend
FeatureProbes = importlib.import_module("features").FeatureProbes
generate_catalog = importlib.import_module("inventory").generate_catalog

arguments = argparse.ArgumentParser(add_help=False)
arguments.add_argument("--base", type=Path)
options, remaining = arguments.parse_known_args()


class ModuleOptions(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory(prefix="plumed-module-options-")
        self.root = Path(self.scratch.name).resolve()
        self.addCleanup(self.scratch.cleanup)

    def policy(self, text=""):
        (self.root / "plugin.yml").write_text(text or "conda: {}\n")
        return read_policy(self.root)

    def test_defaults_and_missing_file(self):
        policy = read_policy(self.root)
        self.assertFalse(policy["dummy_package"])
        self.assertEqual(policy["test_hardware"], "none")
        self.assertFalse(policy["test_force_requirements"])
        self.assertEqual(policy["required_features"], [])

    def test_kernel_identity_changes_with_module_selection(self):
        tree = self.root / "identity"
        config = tree / "src/config"
        config.mkdir(parents=True)
        tools = tree / "src/maketools"
        tools.mkdir()
        for name in (
            "module-build-info",
            "find_modules.sh",
            "module-selection",
            "module-header",
            "module-header.awk",
        ):
            shutil.copy2(source / "src/maketools" / name, tools / name)
        (tree / "VERSION.txt").write_text("2.11\n")
        (tree / "Makefile.conf").write_text(
            "CXX=c++\nPLUMED_MODULE_ENABLED=\nPLUMED_MODULE_DISABLED=\n"
        )
        for module, kind in (("core", "always"), ("alpha", "default-off")):
            directory = tree / "src" / module
            directory.mkdir()
            (directory / "module.yml").write_text(f"type: {kind}\nuse: []\n")
            (directory / "Source.cpp").write_text("int value;\n")

        def identity():
            subprocess.run(
                ["bash", str(tools / "module-build-info")], cwd=config, check=True
            )
            return (config / "module-build-id").read_text()

        baseline = identity()
        configuration = tree / "Makefile.conf"
        original_config = configuration.read_text()
        configuration.write_text(
            original_config.replace(
                "PLUMED_MODULE_ENABLED=", "PLUMED_MODULE_ENABLED=alpha"
            )
        )
        self.assertNotEqual(identity(), baseline)
        configuration.write_text(original_config)
        self.assertEqual(identity(), baseline)
        declaration = tree / "src/core/module.yml"
        original = declaration.read_text()
        declaration.write_text(
            original
            + "description: Updated prose\nconda: {requirements: {host: [example]}}\n"
        )
        self.assertEqual(identity(), baseline)
        declaration.write_text(original.replace("use: []", "use: [alpha]"))
        self.assertNotEqual(identity(), baseline)
        declaration.write_text(original.replace("type: always", "type: default-on"))
        self.assertNotEqual(identity(), baseline)

    def test_inventory_reads_each_declaration_once(self):
        counts = {}
        original = Path.read_text

        def read(path, *args, **kwargs):
            if path.name in ("module.yml", "plugin.yml"):
                counts[path] = counts.get(path, 0) + 1
            return original(path, *args, **kwargs)

        with mock.patch.object(Path, "read_text", read):
            catalog = generate_catalog(write=False)
        self.assertEqual(len(counts), len(catalog["modules"]) + len(catalog["plugins"]))
        self.assertTrue(all(count == 1 for count in counts.values()))

    def test_complete_unavailable_action_catalogue(self):
        output = self.root / "ModuleMap.inc"
        subprocess.run(
            ["bash", str(source / "src/maketools/makeModuleMap"), str(output)],
            cwd=source / "src/core",
            check=True,
            capture_output=True,
        )
        text = output.read_text()
        for action, module in (
            ("EMMIVOX", "isdb-pytorch"),
            ("DRR", "drr"),
            ("CUDACOORDINATION", "plugins/cudaCoord"),
            ("PYCVINTERFACE", "plugins/pycv"),
            ("PYFUNCTION", "plugins/pycv"),
            ("DISTANCE_VECTORACC", "plugins/openaccPTM"),
        ):
            self.assertIn('"' + action + '", "' + module + '"', text)

    def test_requirements_specs_and_selectors(self):
        policy = self.policy(
            "conda:\n  requirements:\n    host:\n    - hdf5 * nompi_*\n    - spec: libpython\n      selector: py>=315\n    build:\n    - cmake\n    - ninja\n"
        )
        self.assertEqual(
            policy["host"],
            [{"spec": "hdf5 * nompi_*"}, {"spec": "libpython", "selector": "py>=315"}],
        )
        self.assertEqual(policy["build"], [{"spec": "cmake"}, {"spec": "ninja"}])

    def test_malformed_and_unknown_settings(self):
        for text in (
            "conda:\n  host_requirments:\n  - hdf5\n",
            "conda: []",
            'conda:\n  dummy_package: "true"\n',
            "conda:\n  dummy_package: maybe\n",
            "conda:\n  tests:\n    hardware: maybe\n",
            "conda:\n  cxxflags:\n  - spec: -O2\n    selector: linux\n",
            "conda:\n  platforms:\n  - linux\n",
            "conda:\n  required_features:\n  - boost\n",
        ):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.policy(text)

    def test_variant_pairs_and_invalid_combinations(self):
        valid = {
            "variants": {
                "nvhpc_version": ["25.11", "26.1"],
                "openacc_cuda_version": ["13.0", "13.1"],
            },
            "zip_keys": ["nvhpc_version,openacc_cuda_version"],
        }
        policy = self.policy(yaml.safe_dump({"conda": valid}))
        self.assertEqual(policy["variants"]["nvhpc_version"][0]["spec"], "25.11")
        invalid = []
        for key, value in (
            ("openacc_cuda_version", ["13.0"]),
            ("python", []),
            (
                "openacc_cuda_version",
                [
                    {"spec": "13.0", "selector": "linux"},
                    {"spec": "13.1", "selector": "linux"},
                ],
            ),
            ("python", [3.13]),
        ):
            changed = copy.deepcopy(valid)
            changed["variants"][key] = value
            invalid.append(changed)
        invalid.append(dict(valid, zip_keys=["missing,openacc_cuda_version"]))
        for changed in invalid:
            with self.subTest(policy=changed), self.assertRaises(ValueError):
                self.policy(yaml.safe_dump({"conda": changed}))

    def test_variant_template_override(self):
        policy = self.policy(
            "conda:\n  variants:\n    c_stdlib_version:\n    - '2.28'\n    python:\n    - '3.13'\n"
        )
        output = self.root / "conda_build_config.yaml"
        output.with_suffix(".yaml.in").write_text(
            'c_stdlib_version:\n  - "2.17"\nother:\n  - kept\n'
        )
        write_variants(output, [policy])
        text = output.read_text()
        self.assertEqual(text.count("c_stdlib_version:"), 1)
        self.assertNotIn("2.17", text)
        self.assertIn("kept", text)

    def test_shared_matrix_conflict(self):
        first = self.policy("conda:\n  variants:\n    python:\n    - '3.12'\n")
        second = copy.deepcopy(first)
        second["variants"]["python"] = [{"spec": "3.13"}]
        with self.assertRaises(ValueError):
            merge_policies([first, second])

    def test_supported_platforms(self):
        policy = self.policy("conda:\n  platforms:\n  - linux-64\n  - linux-aarch64\n")
        self.assertTrue(supported(policy, "linux-64"))
        self.assertFalse(supported(policy, "osx-arm64"))
        self.assertIn("target_platform == 'linux-64'", skip_selector([policy]))

    def test_backend_restrictions(self):
        policy = self.policy("conda:\n  cmake_args:\n  - -DSOME_OPTION=ON\n")
        validate_backend(policy, "pycv", self.root)
        with self.assertRaises(ValueError):
            validate_backend(policy, "native", self.root)

    def test_environment_arguments_without_shell_evaluation(self):
        policy = self.policy(
            "conda:\n  cppflags:\n  - -I${PREFIX}/include path\n  cxxflags:\n  - -DTEST=1\n"
        )
        env = build_environment(policy, {"PREFIX": str(self.root), "CXXFLAGS": "-O2"})
        self.assertEqual(env["CXXFLAGS"], "-O2 -DTEST=1")
        self.assertIn("include path", env["CPPFLAGS"])
        self.assertEqual(env["PLUMED_CONDA_BUILD"], "1")
        with self.assertRaises(KeyError):
            build_environment(policy, {})

    def test_real_configure_backend_arguments(self):
        policy = self.policy(
            "conda:\n  configure_flags:\n  - --example=${PREFIX}/with space\n  cxxflags:\n  - -DTEST=1\n"
        )
        script = self.root / "configure"
        record = self.root / "record.json"
        script.write_text(
            '#!/usr/bin/env python3\nimport json,os,sys\nfrom pathlib import Path\nPath(os.environ["RECORD"]).write_text(json.dumps([sys.argv[1:],os.environ["CXXFLAGS"]]))\n'
        )
        script.chmod(0o755)
        catalog = {
            "modules": {"tools": {"type": "always", "source_use": [], "conda": policy}},
            "plugins": {},
        }
        self.assertEqual(
            run_build(
                catalog,
                "kernel",
                [str(script)],
                dict(
                    os.environ,
                    PREFIX=str(self.root),
                    RECORD=str(record),
                    CXXFLAGS="-O2",
                ),
            ),
            0,
        )
        self.assertEqual(
            json.loads(record.read_text()),
            [["--example=" + str(self.root / "with space")], "-O2 -DTEST=1"],
        )

    def test_real_make_backend_preserves_existing_flags(self):
        policy = self.policy(
            "conda:\n  cxxflags:\n  - -DTEST=1\n  ldflags:\n  - -Wl,--as-needed\n  make_args:\n  - EXTRA=selected\n"
        )
        catalog = {"modules": {}, "plugins": {"openacc_ptm": {"conda": policy}}}
        (self.root / "Makefile").write_text(
            'include Makefile.conf\nall:\n\t@echo "$(CXXFLAGS)|$(NVLDSHARED)|$(EXTRA)" > recorded\n'
        )
        (self.root / "Makefile.conf").write_text("CXXFLAGS=-O2\nNVLDSHARED=linker\n")
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            # Set target explicitly so this backend unit test runs on macOS too.
            env = dict(os.environ, target_platform="linux-64")
            self.assertEqual(
                run_build(
                    catalog, "openacc_ptm", ["make", "NVLDSHARED=custom-linker"], env
                ),
                0,
            )
        finally:
            os.chdir(previous)
        self.assertEqual(
            (self.root / "recorded").read_text().strip(),
            "-O2 -DTEST=1|custom-linker -Wl,--as-needed|selected",
        )

    def test_required_optional_disabled_probes(self):
        probe = object.__new__(FeatureProbes)
        called = []

        def unavailable(feature):
            called.append(feature)
            return False, [], []

        probe.probe = unavailable
        flags, libraries, features = probe.module_flags(["FFTW"])
        self.assertFalse(features["FFTW"])
        self.assertEqual(libraries, [])
        with self.assertRaisesRegex(ValueError, "Required conda feature FFTW"):
            probe.module_flags(["FFTW"], required=["FFTW"])
        with self.assertRaisesRegex(ValueError, "No native compile/link probe"):
            probe.module_flags([], required=["UNKNOWN"])
        called.clear()
        flags, _, features = probe.module_flags(["FFTW"], disabled=["FFTW"])
        self.assertEqual(called, [])
        self.assertEqual(flags, ["-U__PLUMED_HAS_FFTW"])
        self.assertFalse(features["FFTW"])

    def test_feature_policy_conflict(self):
        with self.assertRaises(ValueError):
            self.policy(
                "conda:\n  required_features:\n  - FFTW\n  disabled_features:\n  - FFTW\n"
            )
        with self.assertRaises(ValueError):
            merge_policies(
                [
                    self.policy("conda:\n  required_features:\n  - FFTW\n"),
                    self.policy("conda:\n  disabled_features:\n  - FFTW\n"),
                ],
                shared_features=True,
            )

    def test_existing_declarations_and_regression_selection(self):
        catalog = generate_catalog(write=False)
        self.assertEqual(
            [
                v["spec"]
                for v in catalog["plugins"]["cuda_coordination"]["conda"]["variants"][
                    "cuda_compiler_version"
                ]
            ],
            ["12.9", "13.0"],
        )
        self.assertEqual(
            len(catalog["plugins"]["pycv"]["conda"]["variants"]["python"]), 5
        )
        tests = regression_selection(catalog)["tests"]
        self.assertEqual(len(tests), 6)
        fourier = next(case for case in tests if case["module"] == "fourier")
        self.assertIn("plumed-ext-symfunc", fourier["packages"])
        self.assertTrue(all(not case["force_requirements"] for case in tests))
        self.assertTrue(
            all(
                not item["conda"]["required_features"]
                for item in catalog["modules"].values()
            )
        )
        for path in [
            *source.glob("src/*/module.yml"),
            *source.glob("plugins/*/plugin.yml"),
        ]:
            if "conda:" in path.read_text():
                self.assertIn("conda/module-options.md", path.read_text(), path)

    @unittest.skipUnless(options.base, "installed kernel supplied by conda integration")
    def test_real_native_extension_policy(self):
        base = options.base.resolve()
        tree = self.root / "source"
        (tree / "extensions").mkdir(parents=True)
        for file in (source / "extensions").glob("*.py"):
            shutil.copy2(file, tree / "extensions" / file.name)
        shutil.copy2(source / "COPYING.LESSER", tree / "COPYING.LESSER")
        (tree / "src/maketools").mkdir(parents=True)
        shutil.copy2(
            source / "src/maketools/module_header.py",
            tree / "src/maketools/module_header.py",
        )
        for plugin in ("pycv", "cudaCoord", "openaccPTM"):
            (tree / "plugins" / plugin / "src").mkdir(parents=True)
        for module, kind in (
            ("tools", "always"),
            ("policy-fixture", "default-off"),
            ("policy-dependency", "default-off"),
        ):
            directory = tree / "src" / module
            directory.mkdir(parents=True)
            (directory / "module.yml").write_text(f"type: {kind}\nuse: []\n")
            if module == "policy-fixture":
                (directory / "module.yml").write_text(
                    f"type: {kind}\nuse: [tools, policy-dependency]\n"
                )
        # USE must remain a native dependency even if ordinary code imports no
        # symbols from B (e.g. a header-only use); linkers may discard unused DSOs.
        (tree / "src/policy-dependency/Unused.cpp").write_text(
            '#include "Value.hpp"\n#ifdef PLUMED_POLICY_TEST\n#error module-specific flag leaked\n#endif\nextern "C" int dependency_test(){return POLICY_VALUE;}\n'
        )
        value_header = tree / "src/policy-dependency/Value.hpp"
        value_header.write_text("#define POLICY_VALUE 7\n")
        for component in ("blas", "lapack"):
            destination = tree / "src" / component
            destination.mkdir()
            shutil.copy2(source / "src" / component / (component + ".h"), destination)
        module = tree / "src/policy-fixture"
        (module / "Check.cpp").write_text(
            '#include "policy-dependency/Value.hpp"\n#ifndef PLUMED_POLICY_TEST\n#error conda flag missing\n#endif\n#ifdef __PLUMED_HAS_FFTW\n#error disabled feature still visible\n#endif\nextern "C" int policy_test(){return POLICY_VALUE+PLUMED_POLICY_TEST;}\n'
        )
        with (module / "module.yml").open("a") as declaration:
            declaration.write(
                yaml.safe_dump(
                    {
                        "conda": {
                            "cxxflags": ["-DPLUMED_POLICY_TEST=1"],
                            "disabled_features": ["FFTW"],
                            "ldflags": [
                                "-Wl,-dead_strip_dylibs"
                                if sys.platform == "darwin"
                                else "-Wl,--as-needed"
                            ],
                        }
                    }
                )
            )
        prefix = self.root / "installed"
        program = base / "bin/plumed"
        # The installed runtime deliberately has no compiler package. This tiny
        # C export/identity fixture uses the test host compiler, retaining the
        # installed ABI CPPFLAGS and testing the module-owned additions.
        compiler_marker = self.root / "host-compiler-used"
        compiler_wrapper = self.root / "host-cxx"
        host_compiler = shlex.split(os.environ.get("CXX", "c++"))
        compiler_wrapper.write_text(
            f"#!{sys.executable}\n"
            "import pathlib,subprocess,sys\n"
            f"pathlib.Path({str(compiler_marker)!r}).touch()\n"
            f"raise SystemExit(subprocess.call({host_compiler!r} + sys.argv[1:]))\n"
        )
        compiler_wrapper.chmod(0o755)
        env = dict(
            os.environ,
            PLUMED_AUTOLOAD="0",
            CXX=str(compiler_wrapper),
            CXXFLAGS="-O2 -std=c++17 -fPIC",
        )
        configuration = subprocess.check_output(
            [str(program), "config", "makefile_conf"], env=env, text=True
        )
        identity = subprocess.check_output(
            [str(program), "--module-build-id"], env=env, text=True
        ).strip()
        soext = next(
            line.split("=", 1)[1]
            for line in configuration.splitlines()
            if line.startswith("SOEXT=")
        )
        kernel = base / "lib" / ("libplumedKernel." + soext)
        kernel_digest = hashlib.sha256(kernel.read_bytes()).digest()
        shutil.copy2(source / "configure", tree / "configure")
        shutil.copy2(source / "Makefile", tree / "Makefile")
        shutil.copy2(source / "extensions/Makefile", tree / "extensions/Makefile")
        env["PATH"] = os.pathsep.join((str(Path(sys.executable).parent), env["PATH"]))
        configure = [
            "bash",
            "configure",
            f"--with-extension-base={base}",
            "--enable-extension-modules=policy-fixture",
            f"--prefix={prefix}",
        ]
        configured = subprocess.run(
            configure,
            cwd=tree,
            env=env,
            text=True,
            capture_output=True,
            check=False,
            timeout=60,
        )
        self.assertEqual(
            configured.returncode, 0, configured.stdout + configured.stderr
        )
        persisted = (tree / "Makefile.conf").read_bytes()
        self.assertTrue(persisted.startswith(configuration.encode()))
        self.assertIn(b"EXTENSION_MODULES=policy-fixture\n", persisted)
        self.assertIn(f"EXTENSION_BUILD_ID={identity}".encode(), persisted)
        for name in ("nonexistent-test-module", "tools"):
            invalid = subprocess.run(
                [
                    "bash",
                    "configure",
                    f"--with-extension-base={base}",
                    f"--enable-extension-modules={name}",
                    f"--prefix={prefix}",
                ],
                cwd=tree,
                env=env,
                text=True,
                capture_output=True,
                check=False,
                timeout=60,
            )
            self.assertNotEqual(invalid.returncode, 0, invalid.stdout + invalid.stderr)
            self.assertEqual((tree / "Makefile.conf").read_bytes(), persisted)
        command = [sys.executable, "extensions/build.py", "install"]
        negative = subprocess.run(
            command,
            check=False,
            cwd=tree,
            env=dict(env, PLUMED_CONDA_BUILD="0"),
            text=True,
            capture_output=True,
            timeout=120,
        )
        self.assertNotEqual(negative.returncode, 0)
        compiler_marker.unlink(missing_ok=True)
        # Makefile.conf assignments override environment variables. Use Make's
        # command-line precedence to select the test host compiler explicitly.
        positive = subprocess.run(
            ["make", "install", f"CXX={env['CXX']}", f"CXXFLAGS={env['CXXFLAGS']}"],
            check=False,
            cwd=tree,
            env=dict(env, PLUMED_CONDA_BUILD="1"),
            text=True,
            capture_output=True,
            timeout=120,
        )
        self.assertEqual(positive.returncode, 0, positive.stdout + positive.stderr)
        self.assertTrue(
            compiler_marker.exists(), "Make ignored the host compiler override"
        )
        self.assertEqual(
            json.loads((tree / "_extensions/modules.json").read_text()),
            ["policy-dependency", "policy-fixture"],
        )
        self.assertEqual(hashlib.sha256(kernel.read_bytes()).digest(), kernel_digest)
        self.assertFalse((prefix / "lib" / kernel.name).exists())
        report = json.loads(
            (prefix / "lib/plumed/plugins/policy-fixture.json").read_text()
        )
        self.assertFalse(report["optional_features"]["FFTW"])
        presence = subprocess.run(
            [
                str(program),
                "--load-module=policy-fixture",
                "modules",
                "--has",
                "policy-fixture",
                "policy-dependency",
                "-q",
            ],
            env=dict(env, PLUMED_PLUGIN_PATH=str(prefix / "lib/plumed/plugins")),
            capture_output=True,
            check=False,
            timeout=30,
        )
        self.assertEqual(presence.returncode, 0, presence.stderr)
        query = subprocess.run(
            [
                str(program),
                "--load-module=policy-fixture",
                "modules",
                "--module",
                "policy-fixture",
                "--has-feature",
                "fftw",
                "-q",
            ],
            check=False,
            env=dict(env, PLUMED_PLUGIN_PATH=str(prefix / "lib/plumed/plugins")),
            capture_output=True,
            timeout=30,
        )
        self.assertEqual(query.returncode, 1, query.stderr)
        library = next((prefix / "lib").glob("libplumed-ext-policy-fixture.*"))

        def value():
            return subprocess.check_output(
                [
                    sys.executable,
                    "-c",
                    "import ctypes,sys; print(ctypes.CDLL(sys.argv[1]).policy_test())",
                    str(library),
                ],
                text=True,
                env=env,
            ).strip()

        self.assertEqual(value(), "8")
        value_header.write_text("#define POLICY_VALUE 41\n")
        for target in ("build", "install"):
            rebuilt = subprocess.run(
                [sys.executable, "extensions/build.py", target],
                cwd=tree,
                env=dict(env, PLUMED_CONDA_BUILD="1"),
                capture_output=True,
                text=True,
                check=False,
                timeout=120,
            )
            self.assertEqual(rebuilt.returncode, 0, rebuilt.stdout + rebuilt.stderr)
        self.assertEqual(
            value(),
            "42",
            "Rebuild used a cached object or stale installed dependency header",
        )


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0], *remaining])
