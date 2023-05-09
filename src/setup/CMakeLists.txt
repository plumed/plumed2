#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "setup")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" ON)
ADDMODULETOKERNEL(${module_name}
Load.cpp
Restart.cpp
Units.cpp
)
ADDMODULENEEDS(${module_name}
	core tools
)
ADDMODULEDEPENDENCIES(${module_name}
	core #tools
)
