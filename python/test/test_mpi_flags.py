"""Communicator annotations without MPI or mpi4py build dependencies."""

import array
import ctypes
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import plumed


@unittest.skipIf(
    sys.platform == "win32", "the test kernels are built with a Unix C++ compiler"
)
class TestMPIFlags(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.kernels = []
        source = Path(__file__).with_name("mpi_flags_kernel.cpp.in")
        compiler = shlex.split(os.environ.get("CXX", "c++"))
        for legacy in (0, 1):
            library = Path(cls.temporary.name) / ("kernel" + str(legacy) + ".so")
            shared = "-dynamiclib" if sys.platform == "darwin" else "-shared"
            subprocess.run(
                compiler
                + [
                    "-std=c++11",
                    shared,
                    "-fPIC",
                    "-DLEGACY_KERNEL=" + str(legacy),
                    "-x",
                    "c++",
                    str(source),
                    "-o",
                    str(library),
                ],
                check=True,
            )
            native = ctypes.CDLL(str(library))
            native.mpi_test_last_flags.restype = ctypes.c_size_t
            cls.kernels.append((library, native))

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def ordinary_values(self):
        values = [
            1,
            1.5,
            "wrong",
            array.array("i", [0]),
            array.array("l", [0]),
            array.array("d", [0.0]),
        ]
        if plumed.HAS_NUMPY:
            import numpy as np

            values += [
                np.zeros(1, dtype=dtype) for dtype in (np.intc, np.int_, np.double)
            ]
        return values

    def test_ordinary_values_and_legacy_kernel(self):
        for legacy, (kernel, native) in enumerate(self.kernels):
            with plumed.Plumed(kernel=kernel) as command:
                for value in self.ordinary_values():
                    with self.subTest(legacy=legacy, value=type(value)):
                        if legacy:
                            command.cmd("setMPIComm", value)
                        else:
                            with self.assertRaisesRegex(
                                TypeError, "pointer to MPI_Comm"
                            ):
                                command.cmd("setMPIComm", value)
                        self.assertTrue(native.mpi_test_last_flags() & 0x20000000)
                command.cmd("setMPIComm", None)
                self.assertFalse(native.mpi_test_last_flags() & 0x20000000)

    def test_recognized_mpi4py_path(self):
        # Simulate the runtime protocol only: no MPI headers/linker dependency.
        class Comm:
            storage = ctypes.c_int(17)

        mpi = types.ModuleType("mpi4py.MPI")
        mpi.Comm = Comm
        mpi._addressof = lambda comm: ctypes.addressof(comm.storage)
        package = types.ModuleType("mpi4py")
        package.MPI = mpi
        with patch.dict(sys.modules, {"mpi4py": package, "mpi4py.MPI": mpi}):
            for kernel, native in self.kernels:
                with plumed.Plumed(kernel=kernel) as command:
                    command.cmd("setMPIComm", Comm())
                    self.assertFalse(native.mpi_test_last_flags() & 0x20000000)
                    self.assertEqual(native.mpi_test_last_flags() & 0xFF0000, 0x10000)

    @unittest.skipUnless(
        os.environ.get("PLUMED_TEST_MPIW_KERNEL"), "requires a fresh MPIW kernel"
    )
    def test_real_kernel_rejects_ordinary_values(self):
        with plumed.Plumed(kernel=os.environ["PLUMED_TEST_MPIW_KERNEL"]) as command:
            for value in self.ordinary_values():
                with self.subTest(value=type(value)):
                    with self.assertRaisesRegex(TypeError, "pointer to MPI_Comm"):
                        command.cmd("setMPIComm", value)


if __name__ == "__main__":
    unittest.main()
