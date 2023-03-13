# The numbers calculated by the following python script are compared with those output from this PLUMED input
#
# MOLINFO STRUCTURE=template.pdb
# t1: TORSION ATOMS=@phi-2
# t2: TORSION ATOMS=@psi-2
# t3: TORSION ATOMS=@phi-3
# t4: TORSION ATOMS=@psi-3
# t5: TORSION ATOMS=@phi-4
# t6: TORSION ATOMS=@psi-4
# t7: TORSION ATOMS=@phi-5
# t8: TORSION ATOMS=@psi-5
# t9: TORSION ATOMS=@phi-6
# t10: TORSION ATOMS=@psi-6
# t11: TORSION ATOMS=@phi-7
# t12: TORSION ATOMS=@psi-7
# t13: TORSION ATOMS=@phi-8
# t14: TORSION ATOMS=@psi-8
# t15: TORSION ATOMS=@phi-9
# t16: TORSION ATOMS=@psi-9
# t17: TORSION ATOMS=@phi-10
# t18: TORSION ATOMS=@psi-10
# t19: TORSION ATOMS=@phi-11
# t20: TORSION ATOMS=@psi-11
# t21: TORSION ATOMS=@phi-12
# t22: TORSION ATOMS=@psi-12
# t23: TORSION ATOMS=@phi-13
# t24: TORSION ATOMS=@psi-13
# t25: TORSION ATOMS=@phi-14
# t26: TORSION ATOMS=@psi-14
# t27: TORSION ATOMS=@phi-15
# t28: TORSION ATOMS=@psi-15
# t29: TORSION ATOMS=@phi-16
# t30: TORSION ATOMS=@psi-16
# t31: TORSION ATOMS=@phi-17
# t32: TORSION ATOMS=@psi-17
#    
# PRINT ARG=t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32 FILE=colvar.ref FMT=%8.4f

import unittest
import numpy as np
import plumed
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

def read_xyz(filename):
   xyz = open(filename)
   n_atoms = int(xyz.readline())
   title, trajectory = xyz.readline(), []
   while True :
      atom_type, coordinates = np.zeros(n_atoms).astype(str), np.zeros([n_atoms,3]) 
      for i in range(0,n_atoms) :
          line = xyz.readline()
          atom,x,y,z = line.split()
          atom_type[i]=atom
          coordinates[i,:]=np.array([x,y,z],dtype=np.float64)
      trajectory.append( coordinates )
      nextline = xyz.readline()
      if( nextline=="" ) : break
      c_atoms = int(nextline)
      if( c_atoms!=n_atoms ) : break 
      title = xyz.readline()
   xyz.close()
   return trajectory

def create_plumed_var( p, name, command ):
   p.cmd("readInputLine", name + ": " + command )
   shape = np.zeros( 1, dtype=np.int_ )
   p.cmd("getDataRank " + name, shape )
   data = np.zeros((1))
   p.cmd("setMemoryForData " + name, data )
   return data

class Test(unittest.TestCase):
  def runtest(self):
    os.system('rm -f bck.*')
    # Output to four decimal places only
    np.set_printoptions(precision=4)
    # Read trajectory
    traj = read_xyz("traj.xyz")
    num_frames = len(traj)
    num_atoms = traj[0].shape[0]
    
    # Create arrays for stuff
    box=np.diag(12.41642*np.ones(3,dtype=np.float64))
    virial=np.zeros((3,3),dtype=np.float64)
    masses=np.ones(num_atoms,dtype=np.float64)
    forces=np.random.rand(num_atoms,3)
    charges=np.zeros(num_atoms,dtype=np.float64)
    
    # Create PLUMED object and read input
    p = plumed.Plumed()

    # not really needed, used to check https://github.com/plumed/plumed2/issues/916
    plumed_version = np.zeros(1, dtype=np.intc)
    p.cmd( "getApiVersion", plumed_version)

    p.cmd("setMDEngine","python")
    p.cmd("setTimestep", 1.)
    p.cmd("setKbT", 1.)
    p.cmd("setNatoms",num_atoms)
    p.cmd("setLogFile","test.log")
    p.cmd("init")
    p.cmd("readInputLine","MOLINFO STRUCTURE=template.pdb")
    t1 = create_plumed_var( p, "t1", "TORSION ATOMS=@phi-2" ) 
    t2 = create_plumed_var( p, "t2", "TORSION ATOMS=@psi-2" )
    t3 = create_plumed_var( p, "t3", "TORSION ATOMS=@phi-3" )
    t4 = create_plumed_var( p, "t4", "TORSION ATOMS=@psi-3" )
    t5 = create_plumed_var( p, "t5", "TORSION ATOMS=@phi-4" )
    t6 = create_plumed_var( p, "t6", "TORSION ATOMS=@psi-4" )
    t7 = create_plumed_var( p, "t7", "TORSION ATOMS=@phi-5" )
    t8 = create_plumed_var( p, "t8", "TORSION ATOMS=@psi-5" )
    t9 = create_plumed_var( p, "t9", "TORSION ATOMS=@phi-6" )
    t10 = create_plumed_var( p, "t10", "TORSION ATOMS=@psi-6" )
    t11 = create_plumed_var( p, "t11", "TORSION ATOMS=@phi-7" )
    t12 = create_plumed_var( p, "t12", "TORSION ATOMS=@psi-7" )
    t13 = create_plumed_var( p, "t13", "TORSION ATOMS=@phi-8" )
    t14 = create_plumed_var( p, "t14", "TORSION ATOMS=@psi-8" )
    t15 = create_plumed_var( p, "t15", "TORSION ATOMS=@phi-9" )
    t16 = create_plumed_var( p, "t16", "TORSION ATOMS=@psi-9" )
    t17 = create_plumed_var( p, "t17", "TORSION ATOMS=@phi-10" )
    t18 = create_plumed_var( p, "t18", "TORSION ATOMS=@psi-10" )
    t19 = create_plumed_var( p, "t19", "TORSION ATOMS=@phi-11" )
    t20 = create_plumed_var( p, "t20", "TORSION ATOMS=@psi-11" )
    t21 = create_plumed_var( p, "t21", "TORSION ATOMS=@phi-12" )
    t22 = create_plumed_var( p, "t22", "TORSION ATOMS=@psi-12" )
    t23 = create_plumed_var( p, "t23", "TORSION ATOMS=@phi-13" )
    t24 = create_plumed_var( p, "t24", "TORSION ATOMS=@psi-13" )
    t25 = create_plumed_var( p, "t25", "TORSION ATOMS=@phi-14" )
    t26 = create_plumed_var( p, "t26", "TORSION ATOMS=@psi-14" )
    t27 = create_plumed_var( p, "t27", "TORSION ATOMS=@phi-15" )
    t28 = create_plumed_var( p, "t28", "TORSION ATOMS=@psi-15" )
    t29 = create_plumed_var( p, "t29", "TORSION ATOMS=@phi-16" )
    t30 = create_plumed_var( p, "t30", "TORSION ATOMS=@psi-16" )
    t31 = create_plumed_var( p, "t31", "TORSION ATOMS=@phi-17" )
    t32 = create_plumed_var( p, "t32", "TORSION ATOMS=@psi-17" )
     
    # Read in the correct answers that were calculated directly using PLUMED
    correct_torsions = np.loadtxt("colvar.ref")
    # Open an output file
    of = open("logfile", "w+")
    
    # Now analyze the trajectory
    for step in range(0,num_frames) :
        of.write("RUNNING ANALYSIS FOR STEP " + str(step) + "\n" )
        p.cmd("setStep",step )
        p.cmd("setBox",box )
        p.cmd("setMasses", masses )
        p.cmd("setCharges", charges )
        p.cmd("setPositions", traj[step])
        p.cmd("setForces", forces )
        p.cmd("setVirial", virial )
        p.cmd("calc")
        bias = np.zeros((1),dtype=np.float64)
        p.cmd("getBias", bias )
        variables = np.array([t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32]).ravel()
        zeros = variables - correct_torsions[step,1:]
        for data in zeros :
            self.assertAlmostEqual(data,0.0,places=4)
    of.close()

  def test(self):
    with cd('test/'):
        self.runtest()

if __name__ == "__main__":
    unittest.main()

