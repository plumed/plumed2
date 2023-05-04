#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "manyrestraints")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" OFF)
ADDMODULETOKERNEL(${module_name}
LWalls.cpp
ManyRestraintsBase.cpp
UWalls.cpp
)
ADDMODULEDEPENDENCIES(${module_name}
core tools vesselbase
)
