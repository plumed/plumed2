import os
import sys
import time

import numpy as np
import plumed

def read_xyz(filename):
   # Based on http://becksteinlab.physics.asu.edu/pages/courses/2013/SimBioNano/03/IntroductiontoPython/p03_instructor.html and change to numpy style.
   # It only reads one frame.
   counter = 0
   xyz = open(filename)
   n_atoms = int(xyz.readline())
   atom_type = np.zeros(n_atoms).astype(str)
   coordinates = np.zeros([n_atoms,3])
   title = xyz.readline()
   for line in xyz:
       if ( (counter+1) > n_atoms ):
                raise ValueError("File says %d atoms but the number of atom lines is greater." % (n_atoms))
       atom,x,y,z = line.split()
       atom_type[counter]=atom
       coordinates[counter,:]=np.array([x,y,z],dtype=np.float64)
       counter += 1
   xyz.close()
   return atom_type, coordinates

# Read XYZ
# filename= os.path.join(
#     os.path.dirname(os.path.abspath(__file__)),
#     os.pardir,
#     os.pardir,
#     "regtest", "crystallization", "rt-q6", "64.xyz")
atom_type, pos = read_xyz("methane.xyz")
num_atoms=pos.shape[0]
step=1
box=np.diag(12.41642*np.ones(3,dtype=np.float64))
virial=np.zeros((3,3),dtype=np.float64)
masses=np.ones(num_atoms,dtype=np.float64)
forces=np.random.rand(num_atoms,3)
charges=np.zeros(num_atoms,dtype=np.float64)

p = plumed.Plumed(8)
p.cmd("setMDEngine","python")
p.cmd("setTimestep", 1.) # Not used but must be defined
p.cmd("setKbT", 1.) # Not used but must be defined
p.cmd("setNatoms",num_atoms)
p.cmd("setPlumedDat","plumed.dat") # Empty, will use the 'action' command
# TODO: write to memory, or disable completely logging
p.cmd("setLogFile","test.log")
# Init
p.cmd("init")
#p.cmd("createAction","d1: DISTANCE ATOMS=1,2")
#p.cmd("createAction","PRINT ARG=d1 FILE=colv")

# Now do calculation
p.cmd("setStep",step )
p.cmd("setBox",box )
p.cmd("setMasses", masses )
p.cmd("setCharges", charges )
p.cmd("setPositions", pos )
p.cmd("setForces", forces )
p.cmd("setVirial", virial )
p.cmd("calc")
bias = np.zeros((1),dtype=np.float64)
p.cmd("getBias", bias )
print("HELLO FINAL BIAS",bias )
#dd = p.grab("d1")
#print("DISTANCE ", dd )



