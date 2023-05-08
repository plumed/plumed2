#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "generic")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" ON)
ADDMODULETOKERNEL(${module_name}
Debug.cpp
DumpAtoms.cpp
DumpDerivatives.cpp
DumpForces.cpp
DumpMassCharge.cpp
DumpProjections.cpp
EffectiveEnergyDrift.cpp
EndPlumed.cpp
FitToTemplate.cpp
Flush.cpp
Group.cpp
Include.cpp
MolInfo.cpp
Plumed.cpp
Print.cpp
RandomExchanges.cpp
Read.cpp
ResetCell.cpp
Time.cpp
UpdateIf.cpp
WholeMolecules.cpp
WrapAround.cpp
)
ADDMODULENEEDS(${module_name}
	core tools xdrfile
)
ADDMODULEDEPENDENCIES(${module_name}
	core tools xdrfile
)
