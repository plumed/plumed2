# The numbers calculated by the following python script are compared with those output from this PLUMED input
#
# MOLINFO STRUCTURE=template.pdb
# t: TORSION ... 
#      ATOMS1=@phi-2
#      ATOMS2=@psi-2
#      ATOMS3=@phi-3
#      ATOMS4=@psi-3
#      ATOMS5=@phi-4
#      ATOMS6=@psi-4
#      ATOMS7=@phi-5
#      ATOMS8=@psi-5
#      ATOMS9=@phi-6
#      ATOMS10=@psi-6
#      ATOMS11=@phi-7
#      ATOMS12=@psi-7
#      ATOMS13=@phi-8
#      ATOMS14=@psi-8
#      ATOMS15=@phi-9
#      ATOMS16=@psi-9
#      ATOMS17=@phi-10
#      ATOMS18=@psi-10
#      ATOMS19=@phi-11
#      ATOMS20=@psi-11
#      ATOMS21=@phi-12
#      ATOMS22=@psi-12
#      ATOMS23=@phi-13
#      ATOMS24=@psi-13
#      ATOMS25=@phi-14
#      ATOMS26=@psi-14
#      ATOMS27=@phi-15
#      ATOMS28=@psi-15
#      ATOMS29=@phi-16
#      ATOMS30=@psi-16
#      ATOMS31=@phi-17
#      ATOMS32=@psi-17
#    

import numpy as np
import plumed

def read_xyz(filename):
   xyz = open(filename)
   n_atoms = int(xyz.readline())
   xyz.readline() # skip line
   trajectory = []
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
      xyz.readline() # skip line
   xyz.close()
   return trajectory

def create_plumed_var( p, name, command ):
   p.cmd("readInputLine", name + ": " + command )
   rank = np.zeros( 1, dtype=np.int_ )
   p.cmd("getDataRank " + name, rank )
   shape = np.zeros( rank, dtype=np.int_ )
   p.cmd("getDataShape " + name, shape )
   data = np.zeros( shape )
   p.cmd("setMemoryForData " + name, data )
   return data

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
p.cmd("setMDEngine","python")
p.cmd("setTimestep", 1.)
p.cmd("setKbT", 1.)
p.cmd("setNatoms",num_atoms)
p.cmd("setLogFile","test.log")
p.cmd("init")
p.cmd("readInputLine","MOLINFO STRUCTURE=template.pdb")
p.cmd("readInputLine","tpsi: TORSION ATOMS1=@psi-2 ATOMS2=@psi-3 ATOMS3=@psi-4 ATOMS4=@psi-5 ATOMS5=@psi-6 ATOMS6=@psi-7 ATOMS7=@psi-8 ATOMS8=@psi-9 ATOMS9=@psi-10 ATOMS10=@psi-11 ATOMS11=@psi-12 ATOMS12=@psi-13 ATOMS13=@psi-14 ATOMS14=@psi-15 ATOMS15=@psi-16 ATOMS16=@psi-17")
p.cmd("readInputLine","tphi: TORSION ATOMS1=@phi-2 ATOMS2=@phi-3 ATOMS3=@phi-4 ATOMS4=@phi-5 ATOMS5=@phi-6 ATOMS6=@phi-7 ATOMS7=@phi-8 ATOMS8=@phi-9 ATOMS9=@phi-10 ATOMS10=@phi-11 ATOMS11=@phi-12 ATOMS12=@phi-13 ATOMS13=@phi-14 ATOMS14=@phi-15 ATOMS15=@phi-16 ATOMS16=@phi-17")
g = create_plumed_var( p, "g", "KDE ARG=tpsi,tphi GRID_BIN=50,50 GRID_MIN=-pi,-pi GRID_MAX=pi,pi BANDWIDTH=0.5,0.5") 

rg = np.loadtxt("grid.inpt")[:,2]

# Read in the correct answers that were calculated directly using PLUMED
#correct_torsions = np.loadtxt("colvar.ref")
# Open an output file
of = open("logfile", "w+")

total_g = np.zeros( g.shape )

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
    if step>0 : total_g = total_g + g 

zeros = rg - total_g.flatten()
for z in zeros :
    if abs(z)>1E-4 : of.write("MISMATCH BETWEEN VALUE FROM PLUMED AND VALUE FROM PYTHON " +  str(z) + " \n")

of.close()
