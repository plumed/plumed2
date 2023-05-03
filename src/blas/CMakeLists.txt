#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "blas")
#Note that the macros here require this directory added as a subdir of plumed/src
set(module_${module_name} ON CACHE INTERNAL "always active module blas")
ADDMODULETOKERNEL(${module_name}
blas.cpp
)
ADDMODULEDEPENDENCIES(${module_name}

)
