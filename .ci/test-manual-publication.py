#!/usr/bin/env python3
"""Exercise the real publisher with fixture manuals and a Git stub (no network)."""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile


SOURCE = Path(__file__).resolve().parent
MKDOCS = (SOURCE.parent / "new-manual").is_dir()
PDF = b"%PDF-1.4\nfixture manual\n%%EOF\n"


class PublicationTest(unittest.TestCase):
    def publish(self, branch, owner="plumed", unsupported=False, what="doc"):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        (root / ".ci").mkdir()
        for name in ("push", "archive-manual-pdf.py"):
            shutil.copyfile(SOURCE / name, root / ".ci" / name)
        (root / "CHANGES").mkdir()
        (root / "CHANGES" / (branch + ".md")).write_text(
            "\\plumednotmaintained\n" if unsupported else "Supported\n"
        )
        for directory in (
            "user-doc/html",
            "developer-doc/html",
            "developer-doc/coverage",
            "developer-doc/coverage-libs",
        ):
            (root / directory).mkdir(parents=True, exist_ok=True)
            (root / directory / "index.html").write_text(
                "<html>\n<head>\n</head>\n</html>\n"
            )
        (root / "user-doc/manual.pdf").write_bytes(PDF)
        (root / "user-doc/html/index.html").write_text(
            "<html>\n<head>\n</head>\n<body>\n"
            '<a href="../manual.pdf">Download</a>\n'
            '<a href="https://example.org/manual.pdf">Other manual</a>\n'
            "</body>\n</html>\n"
        )
        # Populate both inputs: release branches publish Doxygen; master uses MkDocs.
        # Include a filename with spaces and a nested path, not just index.html.
        nested = root / "user-doc/html/nested/a page.html"
        nested.parent.mkdir()
        nested.write_text(
            "<head>\n</head>\n<A HREF='../../manual.pdf#page=2'>PDF</A>\n"
        )
        shutil.copytree(root / "user-doc/html", root / "new-manual/site")
        (root / "bin").mkdir()
        git = root / "bin/git"
        git.write_text('#!/bin/sh\nif test "$1" = log; then echo fixture; fi\nexit 0\n')
        git.chmod(0o755)
        env = dict(
            os.environ,
            GIT_TOKEN="fixture-not-a-secret",
            GITHUB_ACTIONS="true",
            GITHUB_REF="refs/heads/" + branch,
            GITHUB_BASE_REF="",
            TRAVIS="",
            GITHUB_REPOSITORY=owner + "/plumed2",
            PATH=str(root / "bin") + os.pathsep + os.environ["PATH"],
        )
        subprocess.run(
            ["bash", ".ci/push", what],
            cwd=root,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
        return root, next(path for path in (root / "tmp").iterdir() if path.is_dir())

    def test_supported_release_keeps_pdf(self):
        root, published = self.publish("v2.9")
        if MKDOCS:
            self.assertFalse((published / "user-doc/manual.pdf").exists())
        else:
            self.assertEqual((published / "user-doc/manual.pdf").read_bytes(), PDF)
        self.assertFalse((published / "user-doc/manual.pdf.zip").exists())
        html = (published / "user-doc/html/index.html").read_text()
        self.assertIn('href="../manual.pdf"', html)
        self.assertNotIn('content="noindex"', html)

    def test_nonindexable_manuals_archive_pdf(self):
        for branch, owner, unsupported in (
            ("v2.9", "plumed", True),
            ("master", "plumed", False),
            ("topic", "someone", False),
        ):
            with self.subTest(branch=branch, owner=owner):
                root, published = self.publish(branch, owner, unsupported)
                self.assertFalse((published / "user-doc/manual.pdf").exists())
                if MKDOCS:
                    self.assertFalse((published / "user-doc/manual.pdf.zip").exists())
                else:
                    with ZipFile(published / "user-doc/manual.pdf.zip") as archive:
                        self.assertEqual(archive.namelist(), ["manual.pdf"])
                        self.assertEqual(archive.read("manual.pdf"), PDF)
                html = (published / "user-doc/html/index.html").read_text()
                suffix = "" if MKDOCS else ".zip"
                self.assertIn('href="../manual.pdf' + suffix + '"', html)
                self.assertIn('href="https://example.org/manual.pdf"', html)
                self.assertIn('content="noindex"', html)
                nested = (published / "user-doc/html/nested/a page.html").read_text()
                fragment = "#page=2" if MKDOCS else ""
                self.assertIn(
                    "HREF='../../manual.pdf" + suffix + fragment + "'", nested
                )
                self.assertIn('content="noindex"', nested)
                self.assertEqual((root / "user-doc/manual.pdf").read_bytes(), PDF)

    def test_legacy_archive_helper(self):
        root, published = self.publish("master")
        # The source tree still has its PDF, even when the publisher uses MkDocs.
        command = ["python3", str(SOURCE / "archive-manual-pdf.py"), str(root)]
        subprocess.run(command, check=True)
        self.assertFalse((root / "user-doc/manual.pdf").exists())
        with ZipFile(root / "user-doc/manual.pdf.zip") as archive:
            self.assertEqual(archive.namelist(), ["manual.pdf"])
            self.assertEqual(archive.read("manual.pdf"), PDF)
        self.assertIn(
            'href="../manual.pdf.zip"', (root / "user-doc/html/index.html").read_text()
        )
        before = (root / "user-doc/manual.pdf.zip").read_bytes()
        subprocess.run(command, check=True)
        self.assertEqual((root / "user-doc/manual.pdf.zip").read_bytes(), before)

    def test_coverage_is_unchanged(self):
        root, published = self.publish("master", what="coverage")
        self.assertFalse((published / "user-doc").exists())
        self.assertTrue((published / "coverage/index.html").is_file())
        self.assertEqual((root / "user-doc/manual.pdf").read_bytes(), PDF)


if __name__ == "__main__":
    unittest.main()
