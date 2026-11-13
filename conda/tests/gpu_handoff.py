#!/usr/bin/env python3
"""Check channel notebook/provenance and GPU orchestration without GPU execution."""

import ast
import hashlib
import importlib
import json
import os
import runpy
import subprocess
import sys
import tempfile
import types
import unittest
import zipfile
from pathlib import Path
from unittest import mock

source = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(source / "extensions"))
source_record = importlib.import_module("conda_provenance").source_record


class GpuHandoff(unittest.TestCase):
    def test_notebook_complete_dry_run(self):
        notebook = json.loads((source / "conda/gpu-validation.ipynb").read_text())
        with tempfile.TemporaryDirectory(prefix="plumed-notebook-") as scratch:
            scope = {}
            first = "".join(notebook["cells"][1]["source"]).replace(
                "/content/plumed-gpu-validation", scratch
            )
            cell_file = Path(scratch) / "notebook-cell.py"
            cell_file.write_text(first)
            scope = runpy.run_path(str(cell_file), init_globals=scope)
            commands, downloads = [], []
            revision = "a" * 40

            def download(url, destination, expected_sha256=None):
                destination.parent.mkdir(parents=True, exist_ok=True)
                if url.endswith("releases/latest"):
                    destination.write_text(
                        json.dumps(
                            {
                                "assets": [
                                    {
                                        "name": "Miniforge3-example-Linux-x86_64.sh"
                                        + suffix,
                                        "browser_download_url": "https://example.invalid/installer"
                                        + suffix,
                                    }
                                    for suffix in ("", ".sha256")
                                ]
                            }
                        )
                    )
                elif url.endswith(".sha256"):
                    destination.write_text(
                        hashlib.sha256(b"installer").hexdigest() + " installer"
                    )
                else:
                    self.assertEqual(
                        expected_sha256, hashlib.sha256(b"installer").hexdigest()
                    )
                    destination.write_bytes(b"installer")
                return destination

            def run(command, log_path, *, env=None, cwd=None, check=True):
                command = list(map(str, command))
                commands.append((command, env))
                log_path.parent.mkdir(parents=True, exist_ok=True)
                log_path.write_text("SIMULATED\n")
                if command[0] == "bash" and "-b" in command:
                    target = scope["conda_prefix"] / "bin/conda"
                    target.parent.mkdir(parents=True)
                    target.touch()
                elif "create" in command:
                    prefix = scope["kernel_prefix"]
                    (prefix / "conda-meta").mkdir(parents=True)
                    (prefix / "conda-meta/plumed-test.json").write_text(
                        json.dumps(
                            {
                                "name": "plumed",
                                "version": "2026.11.13",
                                "build": "exact_kernel",
                                "sha256": "b" * 64,
                            }
                        )
                    )
                    (prefix / "lib/plumed").mkdir(parents=True)
                    (prefix / "lib/plumed/conda-source.json").write_text(
                        json.dumps(
                            {
                                "revision": revision,
                                "repository": "plumed/plumed2",
                                "dirty": False,
                            }
                        )
                    )
                elif command[0] == "git" and command[1] == "init":
                    (scope["source_dir"] / ".git").mkdir()
                elif command[0] == "bash" and command[1].endswith("run-gpu-test.sh"):
                    self.assertEqual(
                        env["GPU_KERNEL_SPEC"], "plumed=2026.11.13=exact_kernel"
                    )
                    (scope["run_dir"] / "summary.json").write_text(
                        json.dumps(
                            {
                                "passed": True,
                                "variants": {"future-variant": {"status": "PASS"}},
                            }
                        )
                    )
                return 0

            scope.update(download=download, run_logged=run)
            colab = types.ModuleType("google.colab")
            colab.files = types.SimpleNamespace(download=downloads.append)
            google = types.ModuleType("google")
            google.colab = colab
            with (
                mock.patch("platform.system", return_value="Linux"),
                mock.patch("platform.machine", return_value="x86_64"),
                mock.patch("shutil.which", return_value="nvidia-smi"),
                mock.patch.dict(sys.modules, {"google": google, "google.colab": colab}),
            ):
                for index in (2, 3, 4, 5, 7):
                    cell_file.write_text("".join(notebook["cells"][index]["source"]))
                    scope = runpy.run_path(str(cell_file), init_globals=scope)
            self.assertTrue(scope["summary"]["passed"])
            self.assertTrue(
                any(
                    "fetch" in command and revision in command
                    for command, _ in commands
                )
            )
            self.assertEqual(len(downloads), 1)
            with zipfile.ZipFile(downloads[0]) as archive:
                self.assertIn("tests/summary.json", archive.namelist())
                self.assertIn("setup/package-source.json", archive.namelist())
                self.assertFalse(
                    any(
                        "conda-meta/" in name or name.endswith(".conda")
                        for name in archive.namelist()
                    )
                )

    def test_notebook_uses_published_package_source(self):
        notebook = json.loads((source / "conda/gpu-validation.ipynb").read_text())
        code = ""
        for cell in notebook["cells"]:
            if cell["cell_type"] == "code":
                text = "".join(cell["source"])
                ast.parse(text)
                code += text
        self.assertNotIn("nightly.link", code)
        self.assertNotIn("2.11.0dev5", code)
        self.assertIn('channel = "plumed/label/test"', code)
        self.assertIn('provenance.get("revision")', code)
        self.assertIn("GPU_KERNEL_SPEC=kernel_spec", code)

    def test_source_archive_preserves_embedded_provenance(self):
        with tempfile.TemporaryDirectory(prefix="plumed-provenance-") as scratch:
            root = Path(scratch)
            self.assertIsNone(source_record(root)["revision"])
            (root / "extensions").mkdir()
            expected = {
                "schema": 1,
                "revision": "a" * 40,
                "repository": "plumed/plumed2",
                "dirty": False,
            }
            (root / "extensions/conda-source.json").write_text(json.dumps(expected))
            self.assertEqual(source_record(root), expected)

    def test_installation_numeric_failures_and_exact_kernel_pin(self):
        plan = json.loads(
            subprocess.check_output(
                [sys.executable, str(source / "conda/gpu-plan.py")], text=True
            )
        )
        self.assertTrue(any(item["cuda"] for item in plan))
        self.assertTrue(any(not item["cuda"] for item in plan))
        for failure in ("none", "install", "numeric"):
            with (
                self.subTest(failure=failure),
                tempfile.TemporaryDirectory(prefix="plumed-gpu-handoff-") as scratch,
            ):
                root = Path(scratch)
                tools = root / "bin"
                tools.mkdir()
                (tools / "nvidia-smi").write_text("#!/bin/sh\nexit 0\n")
                conda = tools / "conda"
                conda.write_text("""#!/usr/bin/env python3
import json,os,sys
from pathlib import Path
assert "CONDA_OVERRIDE_CUDA" not in os.environ
args=sys.argv[1:]
assert "plumed=2026.11.13=exact_kernel" in args
assert "--strict-channel-priority" in args
with Path(os.environ["COMMAND_RECORD"]).open("a") as f: f.write(json.dumps(args)+"\\n")
if os.environ["FAILURE_MODE"]=="install" and any(os.environ["FAIL_BUILD"] in arg for arg in args): sys.exit(2)
prefix=Path(args[args.index("-p")+1]);(prefix/"conda-meta").mkdir(parents=True);(prefix/"bin").mkdir()
(prefix/"conda-meta/plumed-test.json").write_text(json.dumps({"name":"plumed","version":"2026.11.13","build":"exact_kernel"}))
program=prefix/"bin/python"
program.write_text('#!/bin/sh\\n[ "$FAILURE_MODE" = numeric ] && exit 7\\necho SIMULATED-HARNESS\\n')
program.chmod(0o755)
""")
                for file in tools.iterdir():
                    file.chmod(0o755)
                (tools / "python").symlink_to(sys.executable)
                record = root / "commands.jsonl"
                environment = dict(
                    os.environ,
                    PATH=str(tools) + os.pathsep + os.environ["PATH"],
                    VERSION="2026.11.13",
                    GPU_CONDA_CHANNEL="plumed/label/test",
                    GPU_KERNEL_SPEC="plumed=2026.11.13=exact_kernel",
                    GPU_TEST_ROOT=str(root / "results"),
                    FAILURE_MODE=failure,
                    FAIL_BUILD=plan[0]["build"],
                    COMMAND_RECORD=str(record),
                    CONDA_OVERRIDE_CUDA="999",
                )
                result = subprocess.run(
                    ["bash", str(source / "conda/run-gpu-test.sh")],
                    env=environment,
                    text=True,
                    capture_output=True,
                    check=False,
                    timeout=30,
                )
                self.assertTrue(
                    (root / "results/summary.json").exists(),
                    result.stdout + result.stderr,
                )
                summary = json.loads((root / "results/summary.json").read_text())
                self.assertEqual(
                    summary["passed"], failure == "none", result.stdout + result.stderr
                )
                self.assertEqual(result.returncode == 0, failure == "none")
                self.assertEqual(len(summary["variants"]), len(plan))
                self.assertEqual(len(record.read_text().splitlines()), len(plan))
                if failure == "install":
                    self.assertEqual(
                        summary["variants"][plan[0]["name"]]["status"],
                        "INSTALLATION FAILED",
                    )
                if failure == "numeric":
                    self.assertTrue(
                        all(
                            item["status"] == "NUMERICAL TEST FAILED"
                            for item in summary["variants"].values()
                        )
                    )


if __name__ == "__main__":
    unittest.main()
