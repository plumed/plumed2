from __future__ import print_function
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

def test():
    with cd('test/'):
        d=plumed.read_as_pandas("COLVAR")
        print(d,file=open("dataframe","wt"))
        assert filecmp.cmp("dataframe","dataframe.ref")

        d=plumed.read_as_pandas("COLVAR_gzipped.gz")
        print(d,file=open("dataframe","wt"))
        assert filecmp.cmp("dataframe","dataframe.ref")

        i=0
        for d in plumed.read_as_pandas("COLVAR",chunksize=4):
            print(d,file=open("dataframe."+str(i),"wt"))
            assert filecmp.cmp("dataframe."+str(i),"dataframe."+str(i)+".ref")
            i=i+1

        i=0
        for d in plumed.read_as_pandas("COLVAR_gzipped.gz",chunksize=4):
            print(d,file=open("dataframe."+str(i),"wt"))
            assert filecmp.cmp("dataframe."+str(i),"dataframe."+str(i)+".ref")
            i=i+1

if __name__ == "__main__":
    test()

