#!/usr/bin/env python3
"""Check that the importer cannot write outside its extraction directory."""

import importlib.util
import io
from pathlib import Path
import tarfile
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location(
    "thread_mpi_reimport", ROOT / "src/thread_mpi/reimport.py"
)
importer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(importer)


def archive(name, kind=tarfile.REGTYPE):
    data = io.BytesIO()
    with tarfile.open(fileobj=data, mode="w") as tf:
        member = tarfile.TarInfo(name)
        member.type = kind
        if member.isfile():
            member.size = 5
            tf.addfile(member, io.BytesIO(b"hello"))
        else:
            member.linkname = "../outside"
            tf.addfile(member)
    return data.getvalue()


class ExtractionTests(unittest.TestCase):
    def test_regular_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            importer._extract_archive(archive("src/include/example.h"), root)
            self.assertEqual((root / "src/include/example.h").read_bytes(), b"hello")

    def test_unsafe_members(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "extract"
            root.mkdir()
            for name, kind in [
                ("../outside", tarfile.REGTYPE),
                (str(Path(tmp) / "outside"), tarfile.REGTYPE),
                ("src/../../outside", tarfile.REGTYPE),
                ("src/link", tarfile.SYMTYPE),
                ("src/link", tarfile.LNKTYPE),
                ("src/device", tarfile.CHRTYPE),
                ("src/pipe", tarfile.FIFOTYPE),
            ]:
                with self.subTest(name=name, kind=kind):
                    with self.assertRaises(ValueError):
                        importer._extract_archive(archive(name, kind), root)
                    self.assertEqual(list(root.rglob("*")), [])
                    self.assertFalse((Path(tmp) / "outside").exists())

    def test_existing_directory_symlink(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "extract"
            outside = Path(tmp) / "outside"
            root.mkdir()
            outside.mkdir()
            (root / "include").symlink_to(outside, target_is_directory=True)
            with self.assertRaises(ValueError):
                importer._extract_archive(archive("include/example.h"), root)
            self.assertEqual(list(outside.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
