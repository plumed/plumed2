#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "mapping")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" ON)
ADDMODULETOKERNEL(${module_name}
AdaptivePath.cpp
Mapping.cpp
PathBase.cpp
Path.cpp
PathReparameterization.cpp
PathTools.cpp
PCAVars.cpp
PropertyMap.cpp
SpathVessel.cpp
TrigonometricPathVessel.cpp
ZpathVessel.cpp
)
ADDMODULENEEDS(${module_name}
	core tools cltools reference vesselbase
)
ADDMODULEDEPENDENCIES(${module_name}
	cltools #tools core
	reference #tools core
	vesselbase #tools core
)
