#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "molfile")
#Note that the macros here require this directory added as a subdir of plumed/src
message (FATAL_ERROR "molfile has a non standard Makefile (more than 4 lines) you need to modify the CMakeLists.txt!")
option(module_${module_name} "activate module molfile" ON)
ADDMODULETOKERNEL(${module_name}
crdplugin.cpp
dcdplugin.cpp
gromacsplugin.cpp
pdbplugin.cpp
)
ADDMODULEDEPENDENCIES(${module_name}

)
