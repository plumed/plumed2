# Values passed to and returned by cmd whose C types depend on the platform
# (e.g. long is 32 bit on Windows, 64 bit on Linux and macOS)
import array
import os
import tempfile
import unittest

import plumed

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def setUpModule():
    global LOGDIR
    LOGDIR = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)


def tearDownModule():
    LOGDIR.cleanup()


def distance_kernel():
    p = plumed.Plumed()
    p.cmd("setNatoms", 4)
    p.cmd("setLogFile", os.path.join(LOGDIR.name, "test.log"))
    p.cmd("init")
    p.cmd("readInputLine", "d: DISTANCE ATOMS1=1,2 ATOMS2=3,4")
    return p


class TestIntegerArrays(unittest.TestCase):
    @unittest.skipUnless(HAS_NUMPY, "requires numpy")
    def test_numpy_integers(self):
        # np.int_ is 64 bit also on Windows with numpy>=2, where long is 32 bit
        for dtype in (np.int_, np.int64):
            with self.subTest(dtype=np.dtype(dtype).name), distance_kernel() as p:
                rank = np.zeros(1, dtype=dtype)
                p.cmd("getDataRank d", rank)
                self.assertEqual(rank[0], 1)
                shape = np.zeros(1, dtype=dtype)
                p.cmd("getDataShape d", shape)
                self.assertEqual(shape[0], 2)
        with distance_kernel() as p, self.assertRaises(ValueError):
            p.cmd("getDataRank d", np.zeros(1, dtype=np.int16))

    def test_array_long(self):
        with distance_kernel() as p:
            rank = array.array("l", [0])
            p.cmd("getDataRank d", rank)
            self.assertEqual(rank[0], 1)
            shape = array.array("l", [0])
            p.cmd("getDataShape d", shape)
            self.assertEqual(shape[0], 2)


class TestExceptions(unittest.TestCase):
    def test_filesystem_error_paths(self):
        # the kernel passes the paths as native strings without a terminator
        # (wchar_t on Windows): they must be decoded with their length
        with plumed.Plumed() as p, self.assertRaises(OSError) as raised:
            p.cmd(
                "throw",
                "std::filesystem::filesystem_error std::generic_category 2 "
                "first_path.dat second_path_longer.dat",
            )
        self.assertEqual(raised.exception.errno, 2)
        self.assertEqual(raised.exception.filename, "first_path.dat")
        self.assertEqual(raised.exception.filename2, "second_path_longer.dat")

    def test_filesystem_error_without_paths(self):
        with plumed.Plumed() as p, self.assertRaises(OSError) as raised:
            p.cmd("throw", "std::filesystem::filesystem_error std::generic_category 2")
        self.assertIsNone(raised.exception.filename)
        self.assertIsNone(raised.exception.filename2)


if __name__ == "__main__":
    unittest.main()
