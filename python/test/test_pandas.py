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

def test1():
    with cd('test/'):
        d=plumed.read_as_pandas("COLVAR",enable_constants='columns')
        print(d,file=open("dataframe","wt"))
        assert filecmp.cmp("dataframe","dataframe.ref")

def test2():
    with cd('test/'):
        d=plumed.read_as_pandas("COLVAR_gzipped.gz",enable_constants='columns')
        print(d,file=open("dataframe","wt"))
        assert filecmp.cmp("dataframe","dataframe.ref")

def test3():
    with cd('test/'):
        i=0
        for d in plumed.read_as_pandas("COLVAR",chunksize=4,enable_constants='columns'):
            print(d,file=open("dataframe."+str(i),"wt"))
            assert filecmp.cmp("dataframe."+str(i),"dataframe."+str(i)+".ref")
            i=i+1

def test4():
    with cd('test/'):
        i=0
        for d in plumed.read_as_pandas("COLVAR_gzipped.gz",chunksize=4,enable_constants='columns'):
            print(d,file=open("dataframe."+str(i),"wt"))
            assert filecmp.cmp("dataframe."+str(i),"dataframe."+str(i)+".ref")
            i=i+1

def test5():
    with cd('test/'):
        d=plumed.read_as_pandas("COLVAR",enable_constants='metadata')
        print(d,file=open("dataframe_noconst","wt"))
        assert filecmp.cmp("dataframe_noconst","dataframe_noconst.ref")
        assert d.plumed_constants[0][0]=="a"
        assert d.plumed_constants[0][2]=="pi"

def test6():
    with cd('test/'):
        i=0
        for d in plumed.read_as_pandas("COLVAR",chunksize=4,enable_constants='metadata'):
            print(d,file=open("dataframe_noconst."+str(i),"wt"))
            assert filecmp.cmp("dataframe_noconst."+str(i),"dataframe_noconst."+str(i)+".ref")
            assert d.plumed_constants[0][0]=="a"
            assert d.plumed_constants[0][2]=="pi"
            i=i+1

def test7():
    with cd('test/'):
        d=plumed.read_as_pandas("COLVAR")
        plumed.write_pandas(d,"COLVAR_write1")
        assert filecmp.cmp("COLVAR_write1","COLVAR_write.ref")
        d=plumed.read_as_pandas("COLVAR",index_col='time')
        plumed.write_pandas(d,"COLVAR_write2")
        assert filecmp.cmp("COLVAR_write2","COLVAR_write.ref")
        d=plumed.read_as_pandas("COLVAR",index_col=('time','psi'))
        try:
            plumed.write_pandas(d,"COLVAR_write3")
            assert False
        except TypeError:
            pass

if __name__ == "__main__":
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()

