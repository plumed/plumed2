#!/usr/bin/env python3
"""Ensure the CodeQL exceptions cannot hide other diagnostics."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path(__file__).with_name("filter-sarif.py")
spec = importlib.util.spec_from_file_location("codeql_filter", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def result(rule, path, line=1):
    return {
        "ruleId": rule,
        "locations": [
            {
                "physicalLocation": {
                    "artifactLocation": {"uri": path},
                    "region": {"startLine": line},
                }
            }
        ],
    }


class FilterTests(unittest.TestCase):
    def test_only_imported_style_notes(self):
        self.assertTrue(
            module.excluded(
                result("cpp/commented-out-code", "src/thread_mpi/comm.cpp"), ROOT
            )
        )
        for rule, path in [
            ("cpp/buffer-overflow", "src/thread_mpi/comm.cpp"),
            ("cpp/commented-out-code", "src/core/CLToolMain.cpp"),
            ("py/tarslip", "src/thread_mpi/reimport.py"),
        ]:
            with self.subTest(rule=rule, path=path):
                self.assertFalse(module.excluded(result(rule, path), ROOT))

    def test_only_opaque_api_comparison(self):
        lines = (ROOT / "src/main/main.cpp").read_text().splitlines()
        api = next(
            i for i, line in enumerate(lines, 1) if line.strip() == "if(api>=13) {"
        )
        self.assertTrue(
            module.excluded(
                result("cpp/constant-comparison", "src/main/main.cpp", api), ROOT
            )
        )
        for rule, path, line in [
            ("cpp/constant-comparison", "src/main/main.cpp", api + 1),
            ("cpp/buffer-overflow", "src/main/main.cpp", api),
            ("cpp/constant-comparison", "src/core/CLToolMain.cpp", api),
        ]:
            self.assertFalse(module.excluded(result(rule, path, line), ROOT))

    def test_imported_quality_notes_only(self):
        quality = {
            "cpp/unused-static-variable",
            "cpp/fixme-comment",
            "cpp/large-parameter",
        }
        for rule, path in [
            ("cpp/unused-static-variable", "src/gemmi/include/gemmi/util.hpp"),
            ("cpp/fixme-comment", "src/small_vector/small_vector.h"),
            ("cpp/large-parameter", "src/metatomic/vesin.cpp"),
        ]:
            with self.subTest(rule=rule, path=path):
                self.assertTrue(module.excluded(result(rule, path), ROOT, quality))
        for rule, path in [
            # security queries and unknown rules stay on in imported code
            ("cpp/path-injection", "src/gemmi/include/gemmi/util.hpp"),
            ("cpp/unused-static-variable", "src/gemmi/include/gemmi/util.hpp"),
            # PLUMED code and the import scripts are not imported code
            ("cpp/large-parameter", "src/metatomic/metatomic.cpp"),
            ("cpp/fixme-comment", "src/gemmi/import.sh"),
            ("cpp/large-parameter", "src/core/ActionWithValue.cpp"),
        ]:
            known = set() if rule == "cpp/unused-static-variable" else quality
            with self.subTest(rule=rule, path=path):
                self.assertFalse(module.excluded(result(rule, path), ROOT, known))

    def test_quality_rules_from_metadata(self):
        run = {
            "tool": {
                "driver": {
                    "rules": [{"id": "a", "properties": {"tags": ["security"]}}]
                },
                "extensions": [
                    {"rules": [{"id": "b", "properties": {"tags": ["quality"]}}]},
                    {"rules": [{"id": "c"}]},
                ],
            }
        }
        self.assertEqual(module.quality_rules(run), {"b", "c"})

    def test_only_plumed_h_false_positives(self):
        path = "src/wrapper/Plumed.h"
        lines = (ROOT / path).read_text().splitlines()
        fold = next(
            i
            for i, line in enumerate(lines, 1)
            if line.strip() == "#if ! defined( __cplusplus) /*{*/"
        )
        sizeof = next(
            i
            for i, line in enumerate(lines, 1)
            if line.strip().startswith("return sizeof(mpiCommPointer(val))")
        )
        load = next(
            i
            for i, line in enumerate(lines, 1)
            if line.strip() == "return dlopen(path,mode);"
        )
        for rule, line in [
            ("cpp/commented-out-code", fold),
            ("cpp/sizeof-side-effect", sizeof),
            ("cpp/uncontrolled-process-operation", load),
        ]:
            self.assertTrue(module.excluded(result(rule, path, line), ROOT))
        for rule, line in [
            ("cpp/commented-out-code", fold + 1),
            ("cpp/commented-out-code", sizeof),
            ("cpp/sizeof-side-effect", fold),
            ("cpp/sizeof-side-effect", sizeof + 1),
            ("cpp/buffer-overflow", fold),
            ("cpp/uncontrolled-process-operation", load + 1),
            ("cpp/uncontrolled-process-operation", sizeof),
            ("cpp/path-injection", load),
        ]:
            with self.subTest(rule=rule, line=line):
                self.assertFalse(module.excluded(result(rule, path, line), ROOT))
        other = result("cpp/sizeof-side-effect", "src/core/PlumedMain.cpp", sizeof)
        self.assertFalse(module.excluded(other, ROOT))
        # LOAD and the MPI provider loader keep the rule
        for other_path in ["src/tools/DLLoader.cpp", "src/mpiwrap/loader.cpp"]:
            other = result("cpp/uncontrolled-process-operation", other_path, load)
            self.assertFalse(module.excluded(other, ROOT))

    def test_only_utf8_fopen_path_injection(self):
        path = "src/tools/FileTools.cpp"
        lines = (ROOT / path).read_text().splitlines()
        fopen = next(
            i
            for i, line in enumerate(lines, 1)
            if line.strip() == "return std::fopen(path.c_str(),mode);"
        )
        rename = next(
            i
            for i, line in enumerate(lines, 1)
            if line.strip() == "return std::rename(from.c_str(),to.c_str());"
        )
        self.assertTrue(
            module.excluded(result("cpp/path-injection", path, fopen), ROOT)
        )
        for rule, other_path, line in [
            ("cpp/path-injection", path, fopen + 1),
            ("cpp/path-injection", path, rename),
            ("cpp/buffer-overflow", path, fopen),
            # direct file accesses elsewhere keep the rule
            ("cpp/path-injection", "src/tools/IFile.cpp", fopen),
            ("cpp/path-injection", "src/cltools/Driver.cpp", fopen),
        ]:
            with self.subTest(rule=rule, path=other_path, line=line):
                self.assertFalse(module.excluded(result(rule, other_path, line), ROOT))

    def test_preserves_report_metadata_and_security_results(self):
        kept = result("cpp/buffer-overflow", "src/thread_mpi/comm.cpp")
        report = {
            "version": "2.1.0",
            "runs": [
                {
                    "tool": {"driver": {"name": "CodeQL"}},
                    "results": [
                        result("cpp/commented-out-code", "src/thread_mpi/comm.cpp"),
                        kept,
                    ],
                }
            ],
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cpp.sarif"
            path.write_text(json.dumps(report))
            subprocess.run([sys.executable, str(SCRIPT), tmp], check=True)
            report["runs"][0]["results"] = [kept]
            self.assertEqual(json.loads(path.read_text()), report)


if __name__ == "__main__":
    unittest.main()
