"""Dataframe I/O, including interoperability with native PLUMED HDF5 tables."""

import io
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import plumed

# pandas (and numpy) may be missing in the wheel tests of new Python versions
try:
    import numpy as np
    import pandas as pd
    from pandas.testing import assert_frame_equal
except ModuleNotFoundError as exc:
    if exc.name not in ("numpy", "pandas"):
        raise
    np = pd = None

try:
    import h5py
except ModuleNotFoundError as exc:
    if exc.name != "h5py":
        raise
    h5py = None


@unittest.skipIf(pd is None, "pandas is not installed")
class TextIO(unittest.TestCase):
    def test_borrowed_stream_and_nondefault_index(self):
        source = io.StringIO("#! FIELDS time value\n0 2\n1 3\n")
        df = plumed.read_as_pandas(source, enable_conversion=False)
        self.assertFalse(source.closed)
        df.index = pd.Index([10, 20], name="row")
        target = io.StringIO()
        plumed.write_pandas(df, target)
        self.assertFalse(target.closed)
        self.assertEqual(
            target.getvalue(), "#! FIELDS row time value\n 10 0 2\n 20 1 3\n"
        )

    def test_text_without_h5py(self):
        original = __import__

        def importing(name, *args, **kwargs):
            if name == "h5py":
                raise ImportError("deliberately unavailable")
            return original(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=importing):
            df = plumed.read_as_pandas(
                io.StringIO("#! FIELDS x\n1\n"), enable_conversion=False
            )
            plumed.write_pandas(df, io.StringIO())
            with self.assertRaisesRegex(ImportError, "requires h5py"):
                plumed.write_pandas(df, io.BytesIO(), format="hdf5")

    def test_chunk_context_and_gzip(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "table.gz"
            plumed.write_pandas(pd.DataFrame({"x": range(5)}), path)
            with plumed.read_as_pandas(
                path, chunksize=2, enable_conversion=False
            ) as reader:
                self.assertEqual(len(next(reader)), 2)
            # File ownership is part of the public I/O contract.
            with self.assertRaises(StopIteration):
                next(reader)
            self.assertEqual(
                len(plumed.read_as_pandas(path, enable_conversion=False)), 5
            )


@unittest.skipIf(pd is None, "pandas is not installed")
@unittest.skipIf(h5py is None, "h5py is optional")
class HDF5IO(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "table.h5"
        self.df = pd.DataFrame(
            {
                "time": [0.0, 0.25, 0.5, 0.75, 1.0],
                "count": np.array([2**53 + 1, 2, 3, 4, 5], dtype=np.int64),
                "label": ["alpha", "β", "two words", "", "end"],
            }
        )

    def read(self, path=None, **kwargs):
        return plumed.read_as_pandas(
            self.path if path is None else path, enable_conversion=False, **kwargs
        )

    def test_typed_roundtrip_and_layout(self):
        plumed.write_pandas(self.df, self.path)
        assert_frame_equal(self.read(enable_constants=False), self.df)
        with h5py.File(self.path) as handle:
            self.assertEqual(handle.attrs["plumed_format"], "plumed-table")
            self.assertEqual(handle.attrs["plumed_format_version"], 1)
            group = handle["block_0"]
            self.assertEqual(group.attrs["plumed_columns"], "time count label")
            self.assertEqual(group["count"].dtype, np.dtype("int64"))
            self.assertEqual(group["time"].maxshape, (None,))
            self.assertEqual(group["label"].asstr(encoding="utf-8")[1], "β")

    def test_constants_and_named_index(self):
        df = plumed.read_as_pandas(
            io.StringIO("#! FIELDS time x\n#! SET pi pi\n#! SET unit nm\n0 1\n1 2\n"),
            enable_conversion=False,
            index_col="time",
        )
        plumed.write_pandas(df, self.path)
        result = self.read(index_col="time")
        assert_frame_equal(pd.DataFrame(result), pd.DataFrame(df))
        self.assertEqual(result.plumed_constants, df.plumed_constants)
        columns = self.read(enable_constants="columns")
        self.assertEqual(list(columns["unit"]), ["nm", "nm"])
        output = io.StringIO()
        plumed.write_pandas(result, output)
        self.assertIn("#! SET pi pi\n", output.getvalue())

    def test_typed_constants_roundtrip(self):
        plumed.write_pandas(self.df, self.path)
        with h5py.File(self.path, "r+") as handle:
            block = handle["block_0"]
            block.attrs["plumed_constants"] = "integer real expression numeric_string"
            block.attrs["integer"] = np.int64(2**53 + 1)
            block.attrs["real"] = 1.2345678901234567
            block.attrs["expression"] = "pi"
            block.attrs["numeric_string"] = "3.0"

        class Kernel:
            def cmd(self, command, result):
                expression = command.removeprefix("convert ")
                try:
                    result[0] = np.pi if expression == "pi" else float(expression)
                except ValueError:
                    pass

        with patch("plumed._io.Plumed", Kernel):
            converted = plumed.read_as_pandas(self.path, kernel=Kernel())
        self.assertIsInstance(converted.plumed_constants[0][1], int)
        self.assertEqual(converted.plumed_constants[0][1], 2**53 + 1)
        self.assertEqual(converted.plumed_constants[2][1], np.pi)
        self.assertEqual(converted.plumed_constants[3][1], 3.0)
        converted_path = self.path.with_name("converted.h5")
        plumed.write_pandas(converted, converted_path)
        with h5py.File(converted_path) as handle:
            self.assertEqual(handle["block_0"].attrs["integer"], 2**53 + 1)

        data = self.read()
        self.assertIsInstance(data.plumed_constants[0][1], int)
        self.assertIsInstance(data.plumed_constants[1][1], float)
        output = self.path.with_name("copy.h5")
        plumed.write_pandas(data, output)
        self.assertEqual(self.read(output).plumed_constants, data.plumed_constants)
        with h5py.File(output) as handle:
            self.assertEqual(handle["block_0"].attrs["integer"], 2**53 + 1)
            self.assertEqual(handle["block_0"].attrs["real"], 1.2345678901234567)
        output.write_bytes(b"preserve")
        for constants in ([("large", 2**63)], [("invalid", "\ud800")]):
            data.plumed_constants = constants
            with self.assertRaises((OverflowError, UnicodeError)):
                plumed.write_pandas(data, output)
            self.assertEqual(output.read_bytes(), b"preserve")

    def test_selection(self):
        plumed.write_pandas(self.df, self.path)
        for usecols in (
            [0, 1],
            np.array([0, 1]),
            ["count", "time"],
            lambda name: name != "label",
        ):
            result = self.read(
                usecols=usecols,
                skiprows=[1, 3],
                nrows=2,
                index_col="time",
                enable_constants=False,
            )
            expected = self.df.iloc[[0, 2], :2].set_index("time")
            assert_frame_equal(result, expected)
        result = self.read(skiprows=lambda i: i % 2 == 0, enable_constants=False)
        assert_frame_equal(result, self.df.iloc[[1, 3]].reset_index(drop=True))
        self.assertEqual(len(self.read(skiprows=3)), 2)
        self.assertEqual(len(self.read(nrows=0)), 0)
        for index_col in (np.array(["time", "count"]), pd.Index(["time", "count"])):
            result = self.read(index_col=index_col, enable_constants=False)
            assert_frame_equal(result, self.df.set_index(["time", "count"]))
        assert_frame_equal(self.read(index_col=[], enable_constants=False), self.df)
        original = h5py.Dataset.__getitem__
        sizes = []

        def observe(dataset, key, *args, **kwargs):
            values = original(dataset, key, *args, **kwargs)
            sizes.append(len(values))
            return values

        for count in (0, 1):
            sizes.clear()
            with patch.object(h5py.Dataset, "__getitem__", observe):
                result = self.read(nrows=count, enable_constants=False)
            self.assertEqual(len(result), count)
            self.assertTrue(sizes)
            self.assertLessEqual(max(sizes), count)
        sizes.clear()
        with patch.object(h5py.Dataset, "__getitem__", observe):
            result = self.read(skiprows=3, nrows=1, enable_constants=False)
        assert_frame_equal(result, self.df.iloc[[3]].reset_index(drop=True))
        self.assertLessEqual(max(sizes), 1)

    def test_chunks_and_ownership(self):
        plumed.write_pandas(self.df, self.path)
        with self.read(chunksize=2, enable_constants=False) as reader:
            frames = list(reader)
            self.assertEqual([len(frame) for frame in frames], [2, 2, 1])
            assert_frame_equal(pd.concat(frames), self.df)
        self.assertFalse(reader.file.id.valid)
        with h5py.File(self.path) as file:
            with self.read(file, chunksize=2) as reader:
                next(reader)
            self.assertTrue(file.id.valid)
            with self.assertRaises(StopIteration):
                next(reader)
            self.read(file)
            self.assertTrue(file.id.valid)
        stream = io.BytesIO()
        plumed.write_pandas(self.df, stream, format="hdf5")
        self.assertFalse(stream.closed)
        stream.seek(0)
        assert_frame_equal(
            self.read(stream, format="hdf5", enable_constants=False), self.df
        )
        self.assertFalse(stream.closed)

    def test_empty_and_nan(self):
        df = pd.DataFrame({"x": pd.Series([], dtype=float)})
        plumed.write_pandas(df, self.path)
        assert_frame_equal(self.read(enable_constants=False), df)
        plumed.write_pandas(pd.DataFrame({"x": [np.nan, np.inf, -np.inf]}), self.path)
        self.assertTrue(np.isnan(self.read().x.iloc[0]))
        with h5py.File(self.path, "w") as handle:
            handle.attrs["plumed_format"] = "plumed-table"
            handle.attrs["plumed_format_version"] = 1
        self.assertTrue(self.read().empty)

    def test_multiple_blocks(self):
        plumed.write_pandas(pd.DataFrame({"x": [1, 2]}), self.path)
        with h5py.File(self.path, "r+") as handle:
            group = handle["block_0"]
            group.attrs["plumed_constants"] = "temperature"
            group.attrs["temperature"] = 300.0
            handle.copy("block_0", "block_1")
            handle["block_1"].attrs["temperature"] = 310.0
        with self.assertRaisesRegex(plumed.FormatError, "constants change"):
            self.read()
        result = self.read(enable_constants="columns")
        self.assertEqual(list(result.temperature), [300, 300, 310, 310])
        self.assertEqual(len(self.read(enable_constants=False)), 4)
        selected = self.read(usecols=iter(["x"]), enable_constants=False)
        self.assertEqual(selected.x.tolist(), [1, 2, 1, 2])
        with self.read(chunksize=10) as reader:
            chunks = list(reader)
        self.assertEqual([df.plumed_constants[0][1] for df in chunks], [300, 310])
        with h5py.File(self.path, "r+") as handle:
            handle["block_1"].move("x", "y")
            handle["block_1"].attrs["plumed_columns"] = "y"
        result = self.read(enable_constants=False)
        self.assertEqual(list(result.columns), ["x", "y"])
        self.assertEqual(result.x.notna().sum(), 2)

    def test_block_order_is_numeric(self):
        plumed.write_pandas(pd.DataFrame({"x": [0]}), self.path)
        with h5py.File(self.path, "r+") as handle:
            for i in range(1, 12):
                handle.copy("block_0", "block_" + str(i))
                handle["block_" + str(i)]["x"][0] = i
        self.assertEqual(list(self.read().x), list(range(12)))

    def test_chunk_selection_and_conversion(self):
        plumed.write_pandas(self.df, self.path)
        with self.read(
            chunksize=2, skiprows=[0, 2], nrows=3, enable_constants=False
        ) as reader:
            chunks = list(reader)
        self.assertEqual([len(chunk) for chunk in chunks], [2, 1])
        assert_frame_equal(
            pd.concat(chunks), self.df.iloc[[1, 3, 4]].reset_index(drop=True)
        )
        source = io.StringIO("#! FIELDS x\n#! SET a pi\n1\n")
        df = plumed.read_as_pandas(source, enable_conversion=False)
        plumed.write_pandas(df, self.path)

        class Kernel:
            def cmd(self, command, result):
                if command == "convert pi":
                    result[0] = np.pi

        with patch("plumed._io.Plumed", Kernel):
            converted = plumed.read_as_pandas(self.path, kernel=Kernel())
        self.assertEqual(converted.plumed_constants[0], ("a", np.pi, "pi"))
        self.assertEqual(converted.x.iloc[0], 1)

    def test_changing_schemas_preserve_large_integers(self):
        large = 2**53 + 1
        plumed.write_pandas(pd.DataFrame({"x": [large]}), self.path)
        with h5py.File(self.path, "r+") as handle:
            handle.copy("block_0", "block_1")
            handle["block_1"].move("x", "y")
            handle["block_1"].attrs["plumed_columns"] = "y"
        result = self.read(enable_constants=False)
        self.assertEqual(result.x.dtype, pd.Int64Dtype())
        self.assertEqual(result.x.iloc[0], large)
        self.assertTrue(pd.isna(result.x.iloc[1]))
        with h5py.File(self.path, "r+") as handle:
            del handle["block_1"]
            handle.copy("block_0", "block_1")
            del handle["block_1"]["x"]
            handle["block_1"].create_dataset("x", data=[0.5])
        result = self.read(enable_constants=False)
        self.assertEqual(result.x.tolist(), [large, 0.5])
        self.assertIsInstance(result.x.iloc[0], int)
        result = self.read(enable_constants=False, index_col="x")
        self.assertEqual(result.index.tolist(), [large, 0.5])
        with h5py.File(self.path, "r+") as handle:
            for i in range(2):
                handle["block_" + str(i)].create_dataset("y", data=[i])
                handle["block_" + str(i)].attrs["plumed_columns"] = "x y"
        result = self.read(enable_constants=False, index_col=["x", "y"])
        self.assertEqual(result.index.tolist(), [(large, 0), (0.5, 1)])

    def test_bad_input_and_arguments(self):
        self.path.write_text("text")
        with self.assertRaises(OSError):
            self.read()
        with h5py.File(self.path, "w"):
            pass
        with self.assertRaises(plumed.FormatError):
            self.read()
        plumed.write_pandas(self.df, self.path)
        for kwargs in (
            {"chunksize": 0},
            {"chunksize": 1.5},
            {"nrows": -1},
            {"skiprows": -1},
            {"index_col": True},
            {"usecols": ["missing"]},
            {"format": "bad"},
        ):
            with self.assertRaises(ValueError):
                self.read(**kwargs)
        with h5py.File(self.path, "r+") as handle:
            handle.attrs["plumed_format_version"] = 2
        with self.assertRaisesRegex(plumed.FormatError, "version"):
            self.read()
        with h5py.File(self.path, "r+") as handle:
            handle.attrs["plumed_format_version"] = 1
            handle["block_0"]["count"].resize((2,))
        with self.assertRaisesRegex(plumed.FormatError, "lengths"):
            self.read()

    def test_invalid_output_preserves_destination(self):
        self.path.write_bytes(b"preserve me")
        invalid = [
            pd.DataFrame({"x": [1 + 2j]}),
            pd.DataFrame({"a/b": [1]}),
            pd.DataFrame({"x": [2**63]}, dtype="uint64"),
            pd.DataFrame({"x": [None]}),
            pd.DataFrame([[1, 2]], columns=["x", "x"]),
            pd.DataFrame({"x": [1]}, index=pd.MultiIndex.from_tuples([(0, 1)])),
        ]
        for df in invalid:
            with self.assertRaises((ValueError, TypeError, OverflowError)):
                plumed.write_pandas(df, self.path)
            self.assertEqual(self.path.read_bytes(), b"preserve me")

    def test_explicit_format_and_borrowed_output(self):
        path = self.path.with_suffix(".data")
        plumed.write_pandas(self.df, path, format="hdf5")
        assert_frame_equal(
            self.read(path, format="hdf5", enable_constants=False), self.df
        )
        with h5py.File(self.path, "w") as handle:
            plumed.write_pandas(self.df, handle)
            self.assertTrue(handle.id.valid)
        self.assertEqual(len(self.read()), 5)

    def test_native_interoperability(self):
        executable = shutil.which("plumed")
        if executable is None:
            self.skipTest("Native plumed executable is unavailable")
        probe = subprocess.run(
            [executable, "--no-mpi", "info", "--configuration"],
            capture_output=True,
            text=True,
            check=True,
        )
        if "__PLUMED_HAS_LIBHDF5" not in probe.stdout:
            self.skipTest("Native PLUMED was built without HDF5")
        source = self.path.with_suffix(".dat")
        source.write_text(
            "#! FIELDS time x\n#! SET min_x -pi\n0 1.1234567890123457\n1 2.5\n"
        )
        subprocess.run(
            [
                executable,
                "--no-mpi",
                "h5convert",
                "--input",
                str(source),
                "--output",
                str(self.path),
            ],
            check=True,
            capture_output=True,
        )
        df = self.read()
        self.assertEqual(df.x.iloc[0], 1.1234567890123457)
        self.assertEqual(df.plumed_constants[0][2], "-pi")
        target = self.path.with_name("python.h5")
        plumed.write_pandas(df, target)
        output = self.path.with_name("converted.dat")
        subprocess.run(
            [
                executable,
                "--no-mpi",
                "h5convert",
                "--input",
                str(target),
                "--output",
                str(output),
            ],
            check=True,
            capture_output=True,
        )
        assert_frame_equal(
            pd.DataFrame(plumed.read_as_pandas(output, enable_conversion=False)),
            pd.DataFrame(df),
            check_dtype=False,
        )
        self.assertIn("#! SET min_x -pi", output.read_text())

        # Native string readers need the same HDF5 string type tag as H5Table,
        # while integer columns must retain values beyond float64 precision.
        typed = type(df)(self.df.iloc[[0, 1, 4]].reset_index(drop=True))
        typed.plumed_constants = plumed.Constants(
            [("integer", 2**53 + 1, str(2**53 + 1))]
        )
        plumed.write_pandas(typed, target)
        subprocess.run(
            [
                executable,
                "--no-mpi",
                "h5convert",
                "--input",
                str(target),
                "--output",
                str(output),
            ],
            check=True,
            capture_output=True,
        )
        assert_frame_equal(
            pd.DataFrame(plumed.read_as_pandas(output, enable_conversion=False)),
            pd.DataFrame(typed),
        )
        self.assertIn("#! SET integer " + str(2**53 + 1), output.read_text())

    def test_python_written_table_can_be_restarted_by_plumed(self):
        """A Python-created table can be reopened, trimmed, and appended natively."""
        executable = shutil.which("plumed")
        if executable is None:
            self.skipTest("Native plumed executable is unavailable")
        probe = subprocess.run(
            [executable, "--no-mpi", "info", "--configuration"],
            capture_output=True,
            text=True,
            check=True,
        )
        if "__PLUMED_HAS_LIBHDF5" not in probe.stdout:
            self.skipTest("Native PLUMED was built without HDF5")

        run_dir = self.path.with_name("native-restart")
        run_dir.mkdir()
        table = run_dir / "COLVAR.h5"
        # Seed rows before, at, and after the restart step using only Python.
        initial = pd.DataFrame({"time": [0, 1, 2], "d": [10.0, 20.0, 30.0]})
        plumed.write_pandas(initial, table)
        (run_dir / "plumed.dat").write_text(
            "d: DISTANCE ATOMS=1,2\nPRINT ARG=d FILE=COLVAR.h5\n"
        )
        trajectory = run_dir / "traj.xyz"
        trajectory.write_text(
            "2\n10 10 10\nX 0 0 0\nX 1 0 0\n2\n10 10 10\nX 0 0 0\nX 2 0 0\n"
        )
        completed = subprocess.run(
            [
                executable,
                "--no-mpi",
                "driver",
                "--plumed",
                "plumed.dat",
                "--ixyz",
                str(trajectory),
                "--restart",
                "--initial-step",
                "1",
            ],
            cwd=run_dir,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        result = plumed.read_as_pandas(run_dir / "COLVAR.h5", enable_conversion=False)
        self.assertEqual(result["time"].tolist(), [0, 1, 2])
        np.testing.assert_allclose(result["d"], [10.0, 1.0, 2.0])


if __name__ == "__main__":
    unittest.main()
