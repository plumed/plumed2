# simple test to check setMPIComm and mpi4py
import plumed
import os
import unittest
try:
  from mpi4py import MPI
  HAS_MPI4PY=True
except ImportError:
  HAS_MPI4PY=False

class Test(unittest.TestCase):
  if HAS_MPI4PY:
    def test(self):

      comm = MPI.COMM_WORLD
      p = plumed.Plumed()
      p.cmd("setNatoms",2)
      p.cmd("setLogFile","test.log")
      p.cmd("setMPIComm",comm)
      p.cmd("init")


if __name__ == "__main__":
    unittest.main()
