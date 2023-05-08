#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "membranefusion")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" OFF)
ADDMODULETOKERNEL(${module_name}
FusionPoreExpansionP.cpp
FusionPoreNucleationP.cpp
MemFusionP.cpp
)
ADDMODULENEEDS(${module_name}
	core reference tools colvar
)
ADDMODULEDEPENDENCIES(${module_name}
	core reference tools colvar
)
