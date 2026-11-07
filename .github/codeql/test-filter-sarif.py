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
            i for i, line in enumerate(lines, 1) if line.strip() == "if(api>=11) {"
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
