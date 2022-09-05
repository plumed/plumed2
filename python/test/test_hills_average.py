import unittest
import plumed
import filecmp
import os
from contextlib import contextmanager

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(newdir)
    try:
        yield
    finally:
        os.chdir(prevdir)

class Test(unittest.TestCase):
  def test1(self):
    with cd('test/'):
        plumed.hills_time_average("HILLS_av", "HILLS_av1")
        self.assertTrue(filecmp.cmp("HILLS_av1","HILLS_av1.ref"))

  def test2(self):
    with cd('test/'):
        plumed.hills_time_average("HILLS_av", "HILLS_av2",t0=1.0,t1=1.4)
        self.assertTrue(filecmp.cmp("HILLS_av2","HILLS_av2.ref"))

  def test3(self):
    with cd('test/'):
        plumed.hills_time_average("HILLS_av", "HILLS_av3",frac0=0.2,frac1=0.8)
        self.assertTrue(filecmp.cmp("HILLS_av3","HILLS_av3.ref"))

  def test4(self):
    with cd('test/'):
        df=plumed.read_as_pandas("HILLS_av")
        df1=plumed.hills_time_average(df)
        df["pp"]=1
        plumed.write_pandas(df1,"HILLS_av4")
        self.assertTrue(filecmp.cmp("HILLS_av4","HILLS_av4.ref"))

  def test5(self):
    with cd('test/'):
        df=plumed.read_as_pandas("HILLS_av")
        plumed.hills_time_average(df,inplace=True)
        plumed.write_pandas(df,"HILLS_av5")
        self.assertTrue(filecmp.cmp("HILLS_av5","HILLS_av5.ref"))

if __name__ == "__main__":
    unittest.main()
