#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "funnel")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" OFF)
ADDMODULETOKERNEL(${module_name}
FPS.cpp
Funnel.cpp
)
ADDMODULENEEDS(${module_name}
	core tools colvar bias
)
ADDMODULEDEPENDENCIES(${module_name}
	bias #tools core
	colvar#tools core
)
