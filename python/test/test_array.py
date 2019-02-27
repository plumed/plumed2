# simple test using raw arrays instead of numpy
import plumed
import os
import array

def test():

    p = plumed.Plumed()
    p.cmd("setNatoms",2)
    p.cmd("setLogFile","test.log")
    p.cmd("init")
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2")
    p.cmd("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1")

    box=array.array('d',[10,0,0,0,10,0,0,0,10])
    virial=array.array('d',[0,0,0,0,0,0,0,0,0])
    masses=array.array('d',[1,1])
    charges=array.array('d',[0,0])
    forces=array.array('d',[0,0,0,0,0,0])
    positions=array.array('d',[0,0,0,1,2,3])

    p.cmd("setStep",0)
    p.cmd("setBox",box )
    p.cmd("setMasses", masses )
    p.cmd("setCharges", charges )
    p.cmd("setPositions", positions )
    p.cmd("setForces", forces )
    p.cmd("setVirial", virial )
    p.cmd("calc")

    bias=array.array('d',[0])

    p.cmd("getBias", bias )

    assert (bias[0] - 7.0)**2<1e-8
    assert (forces[0] - 1.0)**2<1e-8
    assert (forces[1] - 2.0)**2<1e-8
    assert (forces[2] - 3.0)**2<1e-8
    assert (forces[3] + 1.0)**2<1e-8
    assert (forces[4] + 2.0)**2<1e-8
    assert (forces[5] + 3.0)**2<1e-8

    

if __name__ == "__main__":
    test()

