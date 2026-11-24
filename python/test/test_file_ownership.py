"""File ownership and interrupt regressions for the maintained Cython interface."""

import gc
import gzip
import importlib.util
import io
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

import plumed

_io = getattr(plumed, "_io", plumed)
_inputbuilder = getattr(plumed, "_inputbuilder", plumed)
_HAS_PANDAS = importlib.util.find_spec("pandas") is not None


@unittest.skipUnless(_HAS_PANDAS, "pandas is not installed")
class TestFileOwnership(unittest.TestCase):
    header = "#! FIELDS time value\n0 1\n1 2\n"

    @contextmanager
    def opened_streams(self):
        """Retain handles so garbage collection cannot hide missing closes."""
        fix_file = _io._fix_file
        opened = []

        def record(file, mode):
            result = fix_file(file, mode)
            if isinstance(result, tuple):
                (stream, owned) = result
            else:
                (stream, owned) = (result, isinstance(file, str) or result is not file)
            if owned:
                opened.append(("owned", stream))
            return result

        try:
            with patch.object(_io, "_fix_file", record):
                yield opened
        finally:
            for _, stream in opened:
                stream.close()

    def test_owned_reads_close_on_success_and_header_error(self):
        for suffix in ("", ".gz"):
            for contents in (self.header, "invalid header\n"):
                with (
                    self.subTest(suffix=suffix, valid=contents == self.header),
                    tempfile.TemporaryDirectory() as directory,
                ):
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
                                plumed.read_as_pandas(filename, enable_conversion=False)
                        self.assertEqual(len(opened), 1)
                        self.assertTrue(opened[0][1].closed)

    def test_owned_reads_close_on_parser_and_postprocessing_errors(self):
        import pandas as pd

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "COLVAR"
            filename.write_text(self.header)
            for target, name in ((pd, "read_csv"), (_io, "Constants")):
                with self.subTest(failure=name), self.opened_streams() as opened:
                    with (
                        patch.object(target, name, side_effect=ValueError("injected")),
                        self.assertRaisesRegex(ValueError, "injected"),
                    ):
                        plumed.read_as_pandas(str(filename), enable_conversion=False)
                    self.assertTrue(opened[0][1].closed)

    def test_owned_writes_close_and_gzip_output_is_complete(self):
        import pandas as pd

        for suffix in ("", ".gz"):
            with (
                self.subTest(suffix=suffix),
                tempfile.TemporaryDirectory() as directory,
            ):
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
                    with (
                        patch.object(_io, "_fix_file", return_value=(stream, owned)),
                        self.assertRaisesRegex(OSError, "injected write failure"),
                    ):
                        plumed.write_pandas(
                            pd.DataFrame({"value": [1]}), "OUTPUT" if owned else stream
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
            with gzip.open(filename, "wt") as stream, self.opened_streams() as opened:
                plumed.write_pandas(pd.DataFrame({"value": [1, 2]}), stream)
                self.assertEqual(opened, [])
                self.assertFalse(stream.closed)
            with gzip.open(filename, "rt") as stream:
                with self.opened_streams() as opened:
                    data = plumed.read_as_pandas(stream, enable_conversion=False)
                    self.assertEqual(opened, [])
                    self.assertFalse(stream.closed)
                self.assertEqual(data.value.tolist(), [1, 2])
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
                with (
                    self.subTest(suffix=suffix, finish=finish),
                    tempfile.TemporaryDirectory() as directory,
                ):
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
                            self.assertEqual([c.value.tolist() for c in reader], [[2]])
                        elif finish == "close":
                            reader.close()
                        else:
                            with (
                                self.assertRaisesRegex(RuntimeError, "stop early"),
                                reader,
                            ):
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
                with (
                    self.subTest(error=type(error).__name__),
                    self.opened_streams() as opened,
                ):
                    reader = plumed.read_as_pandas(
                        str(filename), chunksize=1, enable_conversion=False
                    )
                    with (
                        patch.object(
                            reader.TextFileReader, "__next__", side_effect=error
                        ),
                        self.assertRaises(type(error)),
                    ):
                        next(reader)
                    self.assertTrue(opened[0][1].closed)
            with self.opened_streams() as opened:
                reader = plumed.read_as_pandas(
                    str(filename), chunksize=1, enable_conversion=False
                )
                with (
                    patch.object(
                        _io, "Constants", side_effect=ValueError("conversion")
                    ),
                    self.assertRaisesRegex(ValueError, "conversion"),
                ):
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


class TestRootDiscovery(unittest.TestCase):
    def test_ordinary_failure_uses_legacy_command(self):
        for error in (RuntimeError("cannot create log"), OSError("cannot open log")):
            with self.subTest(error=type(error).__name__):
                with (
                    patch("tempfile.TemporaryDirectory", side_effect=error),
                    patch(
                        "subprocess.check_output", return_value=b"/legacy/root\n"
                    ) as command,
                ):
                    self.assertEqual(
                        _inputbuilder._guessplumedroot("selected-kernel"),
                        "/legacy/root",
                    )
                command.assert_called_once_with(
                    ["plumed", "--no-mpi", "info", "--root"]
                )

    def test_interrupts_propagate_without_fallback(self):
        for error in (KeyboardInterrupt(), SystemExit(7)):
            with self.subTest(error=type(error).__name__):
                with (
                    patch("tempfile.TemporaryDirectory", side_effect=error),
                    patch("subprocess.check_output") as command,
                    self.assertRaises(type(error)) as raised,
                ):
                    _inputbuilder._guessplumedroot("selected-kernel")
                self.assertIs(raised.exception, error)
                command.assert_not_called()


if __name__ == "__main__":
    unittest.main()
