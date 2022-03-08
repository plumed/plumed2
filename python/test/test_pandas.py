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
        d=plumed.read_as_pandas("COLVAR",enable_constants='columns')
        print(d,file=open("dataframe","wt"))
        self.assertTrue(filecmp.cmp("dataframe","dataframe.ref"))

  def test2(self):
    with cd('test/'):
        d=plumed.read_as_pandas("COLVAR_gzipped.gz",enable_constants='columns')
        print(d,file=open("dataframe","wt"))
        self.assertTrue(filecmp.cmp("dataframe","dataframe.ref"))

  def test3(self):
    with cd('test/'):
        i=0
        for d in plumed.read_as_pandas("COLVAR",chunksize=4,enable_constants='columns'):
            print(d,file=open("dataframe."+str(i),"wt"))
            self.assertTrue(filecmp.cmp("dataframe."+str(i),"dataframe."+str(i)+".ref"))
            i=i+1

  def test4(self):
    with cd('test/'):
        i=0
        for d in plumed.read_as_pandas("COLVAR_gzipped.gz",chunksize=4,enable_constants='columns'):
            print(d,file=open("dataframe."+str(i),"wt"))
            self.assertTrue(filecmp.cmp("dataframe."+str(i),"dataframe."+str(i)+".ref"))
            i=i+1

  def test5(self):
    with cd('test/'):
        d=plumed.read_as_pandas("COLVAR",enable_constants='metadata')
        print(d,file=open("dataframe_noconst","wt"))
        self.assertTrue(filecmp.cmp("dataframe_noconst","dataframe_noconst.ref"))
        self.assertEqual(d.plumed_constants[0][0],"a")
        self.assertEqual(d.plumed_constants[0][2],"pi")

  def test6(self):
    with cd('test/'):
        i=0
        for d in plumed.read_as_pandas("COLVAR",chunksize=4,enable_constants='metadata'):
            print(d,file=open("dataframe_noconst."+str(i),"wt"))
            self.assertTrue(filecmp.cmp("dataframe_noconst."+str(i),"dataframe_noconst."+str(i)+".ref"))
            self.assertEqual(d.plumed_constants[0][0],"a")
            self.assertEqual(d.plumed_constants[0][2],"pi")
            i=i+1

  def test7(self):
    with cd('test/'):
        d=plumed.read_as_pandas("COLVAR")
        plumed.write_pandas(d,"COLVAR_write1")
        self.assertTrue(filecmp.cmp("COLVAR_write1","COLVAR_write.ref"))
        d=plumed.read_as_pandas("COLVAR",index_col='time')
        plumed.write_pandas(d,"COLVAR_write2")
        self.assertTrue(filecmp.cmp("COLVAR_write2","COLVAR_write.ref"))
        d=plumed.read_as_pandas("COLVAR",index_col=('time','psi'))
        try:
            plumed.write_pandas(d,"COLVAR_write3")
            self.assertTrue(False)
        except TypeError:
            pass

if __name__ == "__main__":
    unittest.main()
