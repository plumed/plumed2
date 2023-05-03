#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "tools")
#Note that the macros here require this directory added as a subdir of plumed/src
set(module_${module_name} ON CACHE INTERNAL "always active module tools")
ADDMODULETOKERNEL(${module_name}
Angle.cpp
AtomNumber.cpp
BiasRepresentation.cpp
Brent1DRootSearch.cpp
Citations.cpp
Communicator.cpp
ConjugateGradient.cpp
DLLoader.cpp
DynamicList.cpp
ERMSD.cpp
Exception.cpp
FileBase.cpp
ForwardDecl.cpp
Grid.cpp
h36.cpp
HistogramBead.cpp
IFile.cpp
KernelFunctions.cpp
Keywords.cpp
LatticeReduction.cpp
LinkCells.cpp
Log.cpp
LoopUnroller.cpp
Matrix.cpp
MatrixSquareBracketsAccess.cpp
Minimise1DBrent.cpp
MinimiseBase.cpp
MolDataClass.cpp
MultiValue.cpp
NeighborList.cpp
OFile.cpp
OpenMP.cpp
Pbc.cpp
PDB.cpp
PlumedHandle.cpp
Random.cpp
RMSD.cpp
RootFindingBase.cpp
Stopwatch.cpp
Subprocess.cpp
SwitchingFunction.cpp
Tensor.cpp
Tools.cpp
Torsion.cpp
Tree.cpp
TypesafePtr.cpp
Units.cpp
Vector.cpp
)
ADDMODULEDEPENDENCIES(${module_name}
core lapack lepton
)
