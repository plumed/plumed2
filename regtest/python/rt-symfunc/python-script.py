# The numbers calculated by the following python script are compared with those output from this PLUMED input
#
# cmat: CONTACT_MATRIX GROUP=1-64 SWITCH={CUSTOM R_0=4.5 D_MAX=4.5 FUNC=0.5*(cos(pi*x)+1)} COMPONENTS
# beh2: GSYMFUNC_TWOBODY ...
#    WEIGHT=cmat.w VECTORS1=cmat.x VECTORS2=cmat.y VECTORS3=cmat.z
#    FUNCTION1={FUNC=exp(-(x-3)^2) LABEL=g2}
# ...
# beh3: GSYMFUNC_THREEBODY ...
#     WEIGHT=cmat.w VECTORS1=cmat.x VECTORS2=cmat.y VECTORS3=cmat.z
#     FUNCTION1={FUNC=0.25*(cos(pi*sqrt(rjk)/4.5)+1)*exp(-0.1*(rij+rik+rjk))*(1+2*cos(ajik))^2 LABEL=g4}
# ...
# g4m: PRINT ARG=beh2.g2,beh3.g4 FILE=colvar FMT=%8.4f

import numpy as np
import plumed

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
   dpos = name.find(".")
   if dpos < 0 : 
      ll = name
   else :
      ll = name[0:dpos]
   p.cmd("readInputLine", ll + ": " + command )
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
traj = read_xyz("64.xyz")
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
p.cmd("readInputLine","cmat: CONTACT_MATRIX GROUP=1-64 SWITCH={CUSTOM R_0=4.5 D_MAX=4.5 FUNC=0.5*(cos(pi*x)+1)} COMPONENTS")
p.cmd("readInputLine","cmatr: CUSTOM ARG=cmat.x,cmat.y,cmat.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO")
p.cmd("readInputLine","g2_f: CUSTOM ARG=cmatr,cmat.w FUNC=y*exp(-(x-3)^2) PERIODIC=NO")
p.cmd("readInputLine","ones: ONES SIZE=64")
t1 = create_plumed_var( p, "beh2_g2", "MATRIX_VECTOR_PRODUCT ARG=g2_f,ones" ) 
t2 = create_plumed_var( p, "beh3.g4", "GSYMFUNC_THREEBODY WEIGHT=cmat.w ARG=cmat.x,cmat.y,cmat.z FUNCTION1={FUNC=0.25*(cos(pi*sqrt(rjk)/4.5)+1)*exp(-0.1*(rij+rik+rjk))*(1+2*cos(ajik))^2 LABEL=g4}" )
 
# Read in the correct answers that were calculated directly using PLUMED
correct_values = np.loadtxt("colvar.ref")
# Open an output file
of = open("logfile", "w+")

# Now analyze the trajectory
step=0
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
print( "T1:", len(t1), t1 )
print( "T2:", len(t1), t2 )
variables = np.array([t1,t2]).ravel()
zeros = variables - correct_values[1:]
print("ZEROS:", zeros )
for data in zeros :
    if abs(data)>1E-4 : of.write("MISMATCH BETWEEN VALUE FROM PLUMED AND VALUE FROM PYTHON")

of.close()
