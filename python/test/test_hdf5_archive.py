"""The native TAR profile remains readable without extracting its members."""

import io
from pathlib import Path
import tarfile
import tempfile
import shutil
import subprocess
import unittest

import plumed

try:
    import h5py
    import numpy as np
    import pandas as pd
    from pandas.testing import assert_frame_equal
except ImportError:
    h5py = None


@unittest.skipIf(h5py is None, "optional h5py/numpy/pandas are unavailable")
class ArchiveIO(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "table.h5.tar"
        self.df = pd.DataFrame(
            {
                "time": np.arange(8201, dtype=float) / 8,
                "count": np.resize(
                    np.array([2**53 + 1, -(2**63), 2**63 - 1], dtype=np.int64), 8201
                ),
                "label": ["β with spaces"] * 8201,
            }
        )

    def read(self, path=None, **kwargs):
        return plumed.read_as_pandas(
            self.path if path is None else path,
            enable_conversion=False,
            enable_constants=False,
            **kwargs,
        )

    def test_standard_tar_and_typed_roundtrip(self):
        plumed.write_pandas(self.df, self.path)
        assert_frame_equal(self.read(), self.df)
        with tarfile.open(self.path) as archive:
            self.assertEqual(
                archive.getnames(), ["batch_1.h5", "batch_2.h5", "batch_3.h5"]
            )
            for member in archive:
                self.assertEqual((member.offset_data + member.size) % 4096, 0)
                with h5py.File(
                    io.BytesIO(archive.extractfile(member).read()), "r"
                ) as handle:
                    self.assertEqual(handle.attrs["plumed_format_version"], 1)
        renamed = self.path.with_name("renamed.data")
        self.path.rename(renamed)
        assert_frame_equal(self.read(renamed), self.df)

    def test_ordinary_hdf5_and_text_with_archive_suffix(self):
        ordinary = self.path.with_name("ordinary.h5")
        plumed.write_pandas(self.df, ordinary)
        ordinary.rename(self.path)
        assert_frame_equal(self.read(), self.df)
        self.path.write_text("#! FIELDS x\n1\n")
        result = self.read()
        self.assertEqual(list(result.columns), ["x"])
        self.assertEqual(result["x"].tolist(), [1])

    def test_chunks_selection_and_limits_cross_batch_boundaries(self):
        plumed.write_pandas(self.df, self.path)
        with self.read(
            chunksize=3000, skiprows=3, nrows=7001, usecols=["count"]
        ) as reader:
            frames = list(reader)
        self.assertEqual([len(frame) for frame in frames], [3000, 3000, 1001])
        expected = self.df.iloc[3:7004][["count"]].reset_index(drop=True)
        assert_frame_equal(pd.concat(frames), expected)
        with self.read(chunksize=10) as reader:
            self.assertEqual(len(next(reader)), 10)
        self.assertRaises(StopIteration, next, reader)

    def test_incomplete_tail_and_unknown_profile(self):
        plumed.write_pandas(self.df, self.path)
        with tarfile.open(self.path) as archive:
            second = archive.getmembers()[1]
        data = self.path.read_bytes()
        for cut in (second.offset + 7, second.offset + 511, second.offset_data + 17):
            with self.subTest(cut=cut):
                self.path.write_bytes(data[:cut])
                with self.assertWarnsRegex(RuntimeWarning, "incomplete"):
                    frame = self.read()
                assert_frame_equal(frame, self.df.iloc[:4096])
        for torn in (False, True):
            pending = bytearray(data)
            if torn:
                pending[second.offset + 148 : second.offset + 156] = b"broken!!"
            else:
                pending[second.offset : second.offset + 512] = bytes(512)
            self.path.write_bytes(pending)
            if torn:
                with self.assertWarnsRegex(RuntimeWarning, "damaged/unpublished"):
                    frame = self.read()
            else:
                frame = self.read()
            assert_frame_equal(frame, self.df.iloc[:4096])
        bad = bytearray(data)
        bad[265:277] = b"plumed-h5-v9\0"
        bad[148:156] = b"        "
        checksum = sum(bad[:512])
        bad[148:156] = (f"{checksum:06o}\0 ").encode()
        self.path.write_bytes(bad)
        with self.assertRaisesRegex(plumed.FormatError, "version/profile"):
            self.read()

    def test_recovery_record_and_empty_archive(self):
        plumed.write_pandas(self.df, self.path)
        with tarfile.open(self.path) as archive:
            first, second, _ = archive.getmembers()
            image = archive.extractfile(second).read()
        prefix = first.offset_data + first.size
        record = f"plumed-h5-trim-v1 {prefix} 1 {len(image)} 0 4096\n".encode()
        recovery = Path(str(self.path) + ".plumed-trim")
        recovery.write_bytes(record.ljust(512, b"\0") + image)
        self.path.write_bytes(self.path.read_bytes()[: prefix + 16])
        assert_frame_equal(self.read(), self.df.iloc[:8192])
        recovery.unlink()
        plumed.write_pandas(self.df.iloc[:0], self.path)
        self.assertTrue(self.read().empty)

    def test_invalid_conversion_preserves_existing_archive(self):
        plumed.write_pandas(self.df, self.path)
        before = self.path.read_bytes()
        invalid = self.df.copy()
        invalid["label"] = object()
        with self.assertRaises(TypeError):
            plumed.write_pandas(invalid, self.path)
        self.assertEqual(self.path.read_bytes(), before)

    def test_schema_constants_and_integer_precision_across_batches(self):
        frames = [
            pd.DataFrame(
                {"count": np.array([2**53 + 1], dtype="int64"), "label": ["β"]}
            ),
            pd.DataFrame(
                {"count": np.array([-(2**63)], dtype="int64"), "label": ["two words"]}
            ),
            pd.DataFrame({"label": [3.5], "other": ["changed schema"]}),
        ]
        with tarfile.open(self.path, "w", format=tarfile.USTAR_FORMAT) as archive:
            for i, frame in enumerate(frames):
                frame = (
                    plumed.read_as_pandas(
                        io.StringIO("#! FIELDS unused\n0\n"), enable_conversion=False
                    )
                    .iloc[:0]
                    .drop(columns=["unused"])
                    .reindex(frame.index)
                    .assign(**frame)
                )
                frame.plumed_constants = [
                    ("origin", -7 if i < 2 else 3.5, "-7" if i < 2 else "3.5")
                ]
                image = io.BytesIO()
                plumed.write_pandas(frame, image, format="hdf5")
                data = image.getvalue()
                member = tarfile.TarInfo(f"batch_{i + 1}.h5")
                member.uname = "plumed-h5-v1"
                member.gname = f"{0 if i < 2 else 1}:1"
                member.size = (len(data) + 512 + 4095) // 4096 * 4096 - 512
                archive.addfile(
                    member, io.BytesIO(data + bytes(member.size - len(data)))
                )
        result = self.read()
        self.assertEqual(int(result["count"].iloc[0]), 2**53 + 1)
        self.assertEqual(int(result["count"].iloc[1]), -(2**63))
        with plumed.read_as_pandas(
            self.path, chunksize=2, enable_conversion=False
        ) as reader:
            chunks = list(reader)
        self.assertEqual([len(chunk) for chunk in chunks], [2, 1])
        self.assertEqual(chunks[0].plumed_constants[0][1], -7)
        self.assertEqual(chunks[1].plumed_constants[0][1], 3.5)
        with self.assertRaisesRegex(plumed.FormatError, "constants change"):
            plumed.read_as_pandas(self.path, enable_conversion=False)

    @unittest.skipUnless(
        shutil.which("plumed"), "Native plumed executable is unavailable"
    )
    def test_native_export_of_python_archive(self):
        executable = shutil.which("plumed")
        configuration = subprocess.run(
            [executable, "--no-mpi", "info", "--configuration"],
            capture_output=True,
            text=True,
            check=True,
        )
        if "__PLUMED_HAS_LIBHDF5" not in configuration.stdout:
            self.skipTest("Native PLUMED was built without HDF5")
        plumed.write_pandas(self.df, self.path)
        output = self.path.with_suffix(".h5")
        subprocess.run(
            [
                executable,
                "--no-mpi",
                "h5convert",
                "--input",
                str(self.path),
                "--output",
                str(output),
            ],
            check=True,
            capture_output=True,
        )
        assert_frame_equal(self.read(output), self.df)
