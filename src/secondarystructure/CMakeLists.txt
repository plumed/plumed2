#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "secondarystructure")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" ON)
ADDMODULETOKERNEL(${module_name}
AlphaRMSD.cpp
AntibetaRMSD.cpp
ParabetaRMSD.cpp
SecondaryStructureRMSD.cpp
)
ADDMODULENEEDS(${module_name}
	core tools reference vesselbase
)
ADDMODULEDEPENDENCIES(${module_name}
	reference#tools core
	vesselbase #tools core
)
