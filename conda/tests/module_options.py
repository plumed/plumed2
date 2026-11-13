#!/usr/bin/env python3
"""Exercise module packaging declarations, generated policy and build backends."""

import argparse
import copy
import importlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

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
        (self.root / "conda.mk").write_text(text)
        return read_policy(self.root)

    def test_defaults_and_missing_file(self):
        policy = read_policy(self.root)
        self.assertFalse(policy["dummy_package"])
        self.assertEqual(policy["test_hardware"], "none")
        self.assertEqual(policy["test_force_requirements"], "no")
        self.assertEqual(policy["required_features"], [])

    def test_requirements_quoting_selectors_and_replace(self):
        policy = self.policy(
            'CONDA_HOST_REQUIREMENTS = old\nCONDA_HOST_REQUIREMENTS = "hdf5 * nompi_*"\nCONDA_HOST_REQUIREMENTS += libpython  # [py>=315]\nCONDA_BUILD_REQUIREMENTS = cmake ninja\n'
        )
        self.assertEqual(
            policy["host"],
            [{"spec": "hdf5 * nompi_*"}, {"spec": "libpython", "selector": "py>=315"}],
        )
        self.assertEqual(policy["build"], [{"spec": "cmake"}, {"spec": "ninja"}])

    def test_malformed_and_unknown_settings(self):
        for text in (
            "CONDA_HOST_REQUIRMENTS = hdf5",
            "include other.mk",
            "CONDA_DUMMY_PACKAGE += yes",
            "CONDA_DUMMY_PACKAGE = maybe",
            "CONDA_TEST_HARDWARE = maybe",
            "CONDA_CXXFLAGS = -O2 # [linux]",
            "CONDA_PLATFORMS = linux",
            "CONDA_REQUIRED_FEATURES = boost",
        ):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.policy(text)

    def test_variant_pairs_and_invalid_combinations(self):
        valid = "CONDA_VARIANT_nvhpc_version = 25.11 26.1\nCONDA_VARIANT_openacc_cuda_version = 13.0 13.1\nCONDA_ZIP_KEYS = nvhpc_version,openacc_cuda_version\n"
        policy = self.policy(valid)
        self.assertEqual(policy["variants"]["nvhpc_version"][0]["spec"], "25.11")
        for change in (
            valid.replace("13.0 13.1", "13.0"),
            valid.replace(
                "nvhpc_version,openacc_cuda_version", "missing,openacc_cuda_version"
            ),
            "CONDA_VARIANT_python =\n",
            valid.replace("13.0 13.1", "13.0 13.1 # [linux]"),
        ):
            with self.subTest(text=change), self.assertRaises(ValueError):
                self.policy(change)

    def test_variant_template_override(self):
        policy = self.policy(
            "CONDA_VARIANT_c_stdlib_version = 2.28\nCONDA_VARIANT_python = 3.13\n"
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
        first = self.policy("CONDA_VARIANT_python = 3.12\n")
        second = copy.deepcopy(first)
        second["variants"]["python"] = [{"spec": "3.13"}]
        with self.assertRaises(ValueError):
            merge_policies([first, second])

    def test_supported_platforms(self):
        policy = self.policy("CONDA_PLATFORMS = linux-64 linux-aarch64\n")
        self.assertTrue(supported(policy, "linux-64"))
        self.assertFalse(supported(policy, "osx-arm64"))
        self.assertIn("target_platform == 'linux-64'", skip_selector([policy]))

    def test_backend_restrictions(self):
        policy = self.policy("CONDA_CMAKE_ARGS = -DSOME_OPTION=ON\n")
        validate_backend(policy, "pycv", self.root)
        with self.assertRaises(ValueError):
            validate_backend(policy, "native", self.root)

    def test_environment_arguments_without_shell_evaluation(self):
        policy = self.policy(
            'CONDA_CPPFLAGS = "-I${PREFIX}/include path"\nCONDA_CXXFLAGS = -DTEST=1\n'
        )
        env = build_environment(policy, {"PREFIX": str(self.root), "CXXFLAGS": "-O2"})
        self.assertEqual(env["CXXFLAGS"], "-O2 -DTEST=1")
        self.assertIn("include path", env["CPPFLAGS"])
        self.assertEqual(env["PLUMED_CONDA_BUILD"], "1")
        with self.assertRaises(KeyError):
            build_environment(policy, {})

    def test_real_configure_backend_arguments(self):
        policy = self.policy(
            'CONDA_CONFIGURE_FLAGS = "--example=${PREFIX}/with space"\nCONDA_CXXFLAGS = -DTEST=1\n'
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
            "CONDA_CXXFLAGS = -DTEST=1\nCONDA_LDFLAGS = -Wl,--as-needed\nCONDA_MAKE_ARGS = EXTRA=selected\n"
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
                "CONDA_REQUIRED_FEATURES = FFTW\nCONDA_DISABLED_FEATURES = FFTW\n"
            )
        with self.assertRaises(ValueError):
            merge_policies(
                [
                    self.policy("CONDA_REQUIRED_FEATURES = FFTW\n"),
                    self.policy("CONDA_DISABLED_FEATURES = FFTW\n"),
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
            *source.glob("src/*/conda.mk"),
            *source.glob("plugins/*/conda.mk"),
        ]:
            self.assertIn("conda/module-options.md", path.read_text(), path)

    @unittest.skipUnless(options.base, "installed kernel supplied by conda integration")
    def test_real_native_extension_policy(self):
        base = options.base.resolve()
        tree = self.root / "source"
        (tree / "extensions").mkdir(parents=True)
        for file in (source / "extensions").glob("*.py"):
            shutil.copy2(file, tree / "extensions" / file.name)
        shutil.copy2(source / "COPYING.LESSER", tree / "COPYING.LESSER")
        for plugin in ("pycv", "cudaCoord", "openaccPTM"):
            (tree / "plugins" / plugin / "src").mkdir(parents=True)
        for module, kind in (("tools", "always"), ("policy-fixture", "default-off")):
            directory = tree / "src" / module
            directory.mkdir(parents=True)
            (directory / "module.type").write_text(kind + "\n")
            (directory / "Makefile").write_text(
                "USE=tools\n" if module == "policy-fixture" else "USE=\n"
            )
        for component in ("blas", "lapack"):
            destination = tree / "src" / component
            destination.mkdir()
            shutil.copy2(source / "src" / component / (component + ".h"), destination)
        module = tree / "src/policy-fixture"
        (module / "Check.cpp").write_text(
            '#ifndef PLUMED_POLICY_TEST\n#error conda flag missing\n#endif\n#ifdef __PLUMED_HAS_FFTW\n#error disabled feature still visible\n#endif\nextern "C" int policy_test(){return PLUMED_POLICY_TEST;}\n'
        )
        (module / "conda.mk").write_text(
            "CONDA_CXXFLAGS = -DPLUMED_POLICY_TEST=1\nCONDA_DISABLED_FEATURES = FFTW\n"
        )
        prefix = self.root / "installed"
        program = base / "bin/plumed"
        # The installed runtime deliberately has no compiler package. This tiny
        # C export/identity fixture uses the test host compiler, retaining the
        # installed ABI CPPFLAGS and testing the module-owned additions.
        env = dict(
            os.environ,
            PLUMED_AUTOLOAD="0",
            CXX=os.environ.get("CXX", "c++"),
            CXXFLAGS="-O2 -std=c++17 -fPIC",
        )
        configuration = subprocess.check_output(
            [str(program), "config", "makefile_conf"], env=env, text=True
        )
        identity = subprocess.check_output(
            [str(program), "--module-build-id"], env=env, text=True
        ).strip()
        (tree / "Makefile.conf").write_text(
            configuration
            + f"\nEXTENSION_BASE={base}\nEXTENSION_PREFIX={prefix}\nEXTENSION_MODULES=policy-fixture\nEXTENSION_BUILD_ID={identity}\n"
        )
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
        positive = subprocess.run(
            command,
            check=False,
            cwd=tree,
            env=dict(env, PLUMED_CONDA_BUILD="1"),
            text=True,
            capture_output=True,
            timeout=120,
        )
        self.assertEqual(positive.returncode, 0, positive.stdout + positive.stderr)
        report = json.loads(
            (prefix / "lib/plumed/plugins/policy-fixture.json").read_text()
        )
        self.assertFalse(report["optional_features"]["FFTW"])
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


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0], *remaining])
