"""Regression checks for the Python data-file helpers and native conversions."""

from contextlib import contextmanager
import gc
import gzip
import io
import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import plumed


_io = plumed._io

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


@unittest.skipUnless(_HAS_PANDAS, "pandas is not installed")
class TestFileOwnership(unittest.TestCase):
    header = "#! FIELDS time value\n0 1\n1 2\n"

    @contextmanager
    def opened_streams(self):
        """Retain handles so garbage collection cannot hide missing closes."""
        ordinary_open = open
        compressed_open = gzip.open
        opened = []

        def record(opener, kind, *args, **kwargs):
            stream = opener(*args, **kwargs)
            opened.append((kind, stream))
            return stream

        try:
            with (
                patch.object(
                    _io,
                    "open",
                    lambda *a, **k: record(ordinary_open, "plain", *a, **k),
                    create=True,
                ),
                patch.object(
                    gzip,
                    "open",
                    lambda *a, **k: record(compressed_open, "gzip", *a, **k),
                ),
            ):
                yield opened
        finally:
            for _, stream in opened:
                stream.close()

    def test_owned_reads_close_on_success_and_header_error(self):
        for suffix in ("", ".gz"):
            for contents in (self.header, "invalid header\n"):
                with self.subTest(suffix=suffix, valid=contents == self.header):
                    with tempfile.TemporaryDirectory() as directory:
                        filename = str(Path(directory) / ("COLVAR" + suffix))
                        opener = gzip.open if suffix else open
                        with opener(filename, "wt") as stream:
                            stream.write(contents)
                        with self.opened_streams() as opened:
                            if contents == self.header:
                                data = plumed.read_as_pandas(
                                    filename, enable_conversion=False
                                )
                                self.assertEqual(data.value.tolist(), [1, 2])
                            else:
                                with self.assertRaises(plumed.FormatError):
                                    plumed.read_as_pandas(
                                        filename, enable_conversion=False
                                    )
                            self.assertEqual(len(opened), 1)
                            self.assertTrue(opened[0][1].closed)

    def test_owned_reads_close_on_parser_and_postprocessing_errors(self):
        import pandas as pd

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "COLVAR"
            filename.write_text(self.header)
            for target, name in ((pd, "read_csv"), (_io, "Constants")):
                with self.subTest(failure=name):
                    with self.opened_streams() as opened:
                        with patch.object(
                            target, name, side_effect=ValueError("injected")
                        ):
                            with self.assertRaisesRegex(ValueError, "injected"):
                                plumed.read_as_pandas(
                                    str(filename), enable_conversion=False
                                )
                        self.assertTrue(opened[0][1].closed)

    def test_owned_writes_close_and_gzip_output_is_complete(self):
        import pandas as pd

        for suffix in ("", ".gz"):
            with self.subTest(suffix=suffix):
                with tempfile.TemporaryDirectory() as directory:
                    filename = str(Path(directory) / ("OUTPUT" + suffix))
                    with self.opened_streams() as opened:
                        plumed.write_pandas(pd.DataFrame({"value": [1]}), filename)
                        self.assertEqual(len(opened), 1)
                        self.assertTrue(opened[0][1].closed)
                    opener = gzip.open if suffix else open
                    with opener(filename, "rt") as stream:
                        self.assertEqual(stream.read(), "#! FIELDS value\n 1\n")

    def test_write_error_closes_owned_stream_and_preserves_borrowed_stream(self):
        import pandas as pd

        class FailingStream(io.StringIO):
            def write(self, value):
                raise OSError("injected write failure")

        for owned in (False, True):
            with self.subTest(owned=owned):
                stream = FailingStream()
                try:
                    with patch.object(_io, "open", return_value=stream, create=True):
                        with self.assertRaisesRegex(OSError, "injected write failure"):
                            plumed.write_pandas(
                                pd.DataFrame({"value": [1]}),
                                "OUTPUT" if owned else stream,
                            )
                    self.assertEqual(stream.closed, owned)
                finally:
                    stream.close()

    def test_borrowed_streams_and_stdout_remain_open(self):
        import pandas as pd

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "COLVAR"
            filename.write_text(self.header)
            with filename.open() as stream:
                data = plumed.read_as_pandas(stream, enable_conversion=False)
                self.assertFalse(stream.closed)
            self.assertEqual(data.value.tolist(), [1, 2])
            with filename.open("w") as stream:
                plumed.write_pandas(data, stream)
                self.assertFalse(stream.closed)
                stream.write("# still usable\n")
        with io.StringIO() as stream:
            with patch.object(_io.sys, "stdout", stream):
                plumed.write_pandas(pd.DataFrame({"value": [1]}))
            self.assertFalse(stream.closed)
            self.assertEqual(stream.getvalue(), "#! FIELDS value\n 1\n")

    def test_borrowed_read_stream_remains_open_on_error(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "COLVAR"
            filename.write_text("invalid header\n")
            with filename.open() as stream:
                with self.assertRaises(plumed.FormatError):
                    plumed.read_as_pandas(stream, enable_conversion=False)
                self.assertFalse(stream.closed)

    def test_borrowed_gzip_streams_are_used_without_reopening(self):
        import pandas as pd

        with tempfile.TemporaryDirectory() as directory:
            filename = str(Path(directory) / "COLVAR.gz")
            with gzip.open(filename, "wt") as stream:
                with self.opened_streams() as opened:
                    plumed.write_pandas(pd.DataFrame({"value": [1, 2]}), stream)
                    self.assertEqual(opened, [])
                    self.assertFalse(stream.closed)
            with gzip.open(filename, "rt") as stream:
                with self.opened_streams() as opened:
                    data = plumed.read_as_pandas(stream, enable_conversion=False)
                    self.assertEqual(opened, [])
                    self.assertFalse(stream.closed)
                self.assertEqual(data.value.tolist(), [1, 2])
            # Preserve filename-based text decoding for named raw/binary gzip streams.
            for opener in (open, gzip.open):
                with opener(filename, "rb") as stream:
                    with self.opened_streams() as opened:
                        data = plumed.read_as_pandas(stream, enable_conversion=False)
                        self.assertFalse(stream.closed)
                        self.assertEqual(len(opened), 1)
                        self.assertTrue(opened[0][1].closed)
                    self.assertEqual(data.value.tolist(), [1, 2])

    def test_chunked_reads_close_on_exhaustion_explicit_close_and_context_exit(self):
        for suffix in ("", ".gz"):
            for finish in ("exhaust", "close", "context"):
                with self.subTest(suffix=suffix, finish=finish):
                    with tempfile.TemporaryDirectory() as directory:
                        filename = str(Path(directory) / ("COLVAR" + suffix))
                        opener = gzip.open if suffix else open
                        with opener(filename, "wt") as stream:
                            stream.write(self.header)
                        with self.opened_streams() as opened:
                            reader = plumed.read_as_pandas(
                                filename, chunksize=1, enable_conversion=False
                            )
                            self.assertFalse(opened[0][1].closed)
                            self.assertEqual(next(reader).value.tolist(), [1])
                            self.assertFalse(opened[0][1].closed)
                            if finish == "exhaust":
                                self.assertEqual(
                                    [c.value.tolist() for c in reader], [[2]]
                                )
                            elif finish == "close":
                                reader.close()
                            else:
                                with self.assertRaisesRegex(RuntimeError, "stop early"):
                                    with reader:
                                        raise RuntimeError("stop early")
                            self.assertTrue(opened[0][1].closed)
                            reader.close()
                            with self.assertRaises(StopIteration):
                                next(reader)

    def test_chunked_read_errors_close_owned_file(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "COLVAR"
            filename.write_text(self.header)
            for error in (ValueError("parser failure"), KeyboardInterrupt()):
                with self.subTest(error=type(error).__name__):
                    with self.opened_streams() as opened:
                        reader = plumed.read_as_pandas(
                            str(filename), chunksize=1, enable_conversion=False
                        )
                        with patch.object(
                            reader.TextFileReader, "__next__", side_effect=error
                        ):
                            with self.assertRaises(type(error)):
                                next(reader)
                        self.assertTrue(opened[0][1].closed)
            with self.opened_streams() as opened:
                reader = plumed.read_as_pandas(
                    str(filename), chunksize=1, enable_conversion=False
                )
                with patch.object(
                    _io, "Constants", side_effect=ValueError("conversion")
                ):
                    with self.assertRaisesRegex(ValueError, "conversion"):
                        next(reader)
                self.assertTrue(opened[0][1].closed)

    def test_chunked_read_close_preserves_borrowed_file(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "COLVAR"
            filename.write_text(self.header)
            with filename.open() as stream:
                with plumed.read_as_pandas(
                    stream, chunksize=1, enable_conversion=False
                ) as reader:
                    self.assertEqual(next(reader).value.tolist(), [1])
                self.assertFalse(stream.closed)

    def test_abandoned_chunked_reader_releases_owned_file(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "COLVAR"
            filename.write_text(self.header)
            with self.opened_streams() as opened:
                reader = plumed.read_as_pandas(
                    str(filename), chunksize=1, enable_conversion=False
                )
                del reader
                gc.collect()
                self.assertTrue(opened[0][1].closed)


if __name__ == "__main__":
    unittest.main()
