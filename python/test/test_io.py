"""Regression checks for the Python data-file helpers and native conversions."""

import gzip
import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import plumed
from plumed import _io


_HAS_PANDAS = importlib.util.find_spec("pandas") is not None


class TestIO(unittest.TestCase):
    def test_public_helpers_use_python_implementation(self):
        self.assertEqual(Path(_io.__file__).suffix, ".py")
        for name in (
            "FormatError",
            "Constants",
            "read_as_pandas",
            "write_pandas",
            "hills_time_average",
        ):
            self.assertIs(getattr(plumed, name), getattr(_io, name))

    def test_import_does_not_require_pandas(self):
        script = """
import importlib.abc
import sys

class NoPandas(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'pandas' or fullname.startswith('pandas.'):
            raise AssertionError('pandas must be imported lazily')

sys.meta_path.insert(0, NoPandas())
import plumed
assert 'pandas' not in sys.modules
"""
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run(
                [sys.executable, "-c", script],
                cwd=directory,
                check=True,
                capture_output=True,
            )

    def test_converter_borrows_live_kernel(self):
        with plumed.Plumed() as kernel:
            convert = _io._build_convert_function(kernel)
            self.assertIs(convert.kernel, kernel)
            self.assertAlmostEqual(convert("pi"), math.pi)
            self.assertEqual(convert("unrecognized_symbol"), "unrecognized_symbol")
            del convert
            # Releasing the conversion callable must leave its supplied object usable.
            self.assertAlmostEqual(_io._build_convert_function(kernel)("pi"), math.pi)

    def test_converter_unavailable_kernel_warns(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertWarnsRegex(UserWarning, "cannot load PLUMED instance"):
                self.assertIsNone(
                    _io._build_convert_function(str(Path(directory) / "missing-kernel"))
                )

    def test_constants_preserve_original_spelling_and_existing_values(self):
        calls = []

        def convert(value):
            calls.append(value)
            return math.pi

        constants = plumed.Constants(
            [("angle", "pi"), ("already_converted", 2.0, "2")], convert=convert
        )
        self.assertEqual(
            constants, [("angle", math.pi, "pi"), ("already_converted", 2.0, "2")]
        )
        self.assertEqual(calls, ["pi"])
        with self.assertRaisesRegex(ValueError, "2- or 3-plets"):
            plumed.Constants([("invalid",)], convert=convert)

    @unittest.skipUnless(_HAS_PANDAS, "pandas is not installed")
    def test_disabled_conversion_needs_no_kernel(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "COLVAR"
            filename.write_text("#! FIELDS time value\n#! SET angle pi\n0 1\n1 2\n")
            with patch.object(
                _io, "Plumed", side_effect=AssertionError("unexpected kernel load")
            ):
                data = plumed.read_as_pandas(str(filename), enable_conversion=False)
            self.assertEqual(data["value"].tolist(), [1, 2])
            self.assertEqual(data.plumed_constants, [("angle", "pi", "pi")])

    @unittest.skipUnless(_HAS_PANDAS, "pandas is not installed")
    def test_gzip_chunks_and_roundtrip_preserve_constants(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "COLVAR.gz"
            with gzip.open(filename, "wt") as stream:
                stream.write(
                    "#! FIELDS time value\n#! SET angle pi\n0 1\n1 2\n2 3\n3 4\n4 5\n"
                )
            chunks = list(
                plumed.read_as_pandas(
                    str(filename), chunksize=2, enable_conversion=False
                )
            )
            self.assertEqual([len(chunk) for chunk in chunks], [2, 2, 1])
            self.assertEqual(
                [value for chunk in chunks for value in chunk["value"]],
                [1, 2, 3, 4, 5],
            )
            for chunk in chunks:
                self.assertEqual(chunk.plumed_constants, [("angle", "pi", "pi")])
            original = plumed.read_as_pandas(str(filename), enable_conversion=False)
            output = Path(directory) / "ROUNDTRIP.gz"
            plumed.write_pandas(original, str(output))
            reread = plumed.read_as_pandas(str(output), enable_conversion=False)
            self.assertEqual(reread.to_dict("list"), original.to_dict("list"))
            self.assertEqual(reread.plumed_constants, original.plumed_constants)

    @unittest.skipUnless(_HAS_PANDAS, "pandas is not installed")
    def test_hills_time_and_fraction_windows_leave_input_unchanged(self):
        import pandas as pd

        hills = pd.DataFrame({"time": [0, 1, 2, 3, 4], "height": [10.0] * 5})
        average = plumed.hills_time_average(hills)
        self.assertEqual(average.height.tolist(), [10, 7.5, 5, 2.5, 0])
        average = plumed.hills_time_average(hills, t0=1, t1=3)
        self.assertEqual(average.height.tolist(), [10, 10, 5, 0, 0])
        average = plumed.hills_time_average(hills, frac0=0.4, frac1=1)
        self.assertEqual(average.height.tolist(), [10, 10, 10, 5, 0])
        self.assertEqual(hills.height.tolist(), [10] * 5)
        with self.assertRaisesRegex(ValueError, "simultaneously time and fraction"):
            plumed.hills_time_average(hills, t0=1, frac0=0.4)

    @unittest.skipUnless(_HAS_PANDAS, "pandas is not installed")
    def test_hills_inplace_and_file_output(self):
        import pandas as pd

        hills = pd.DataFrame({"time": [0, 1, 2], "height": [6.0, 6.0, 6.0]})
        result = plumed.hills_time_average(hills, inplace=True)
        self.assertIsNone(result)
        self.assertEqual(hills.height.tolist(), [6, 3, 0])
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "HILLS"
            result = plumed.hills_time_average(hills, str(output))
            self.assertIsNone(result)
            reread = plumed.read_as_pandas(str(output), enable_conversion=False)
            self.assertEqual(reread.height.tolist(), [6, 1.5, 0])
            self.assertEqual(hills.height.tolist(), [6, 3, 0])

    @unittest.skipUnless(_HAS_PANDAS, "pandas is not installed")
    def test_invalid_header_uses_public_format_error(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "COLVAR"
            filename.write_text("invalid header\n")
            with self.assertRaises(plumed.FormatError):
                plumed.read_as_pandas(str(filename), enable_conversion=False)


if __name__ == "__main__":
    unittest.main()
