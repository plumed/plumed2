#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "lapack")
#Note that the macros here require this directory added as a subdir of plumed/src
set(module_${module_name} ON CACHE INTERNAL "always active module ${module_name}")
ADDMODULETOKERNEL(${module_name}
    lapack.cpp
)
ADDMODULENEEDS(${module_name}
    blas
)

ADDMODULEDEPENDENCIES(${module_name}
    blas
)

string(TOUPPER ${module_name} externalLib)
if (useExternal_${externalLib} AND ${BLAS_FOUND})
    find_package(${externalLib})
endif (useExternal_${externalLib} AND ${BLAS_FOUND})

if (${externalLib}_FOUND)
    target_compile_definitions(${module_name} PUBLIC __PLUMED_HAS_EXTERNAL_LAPACK=1)
    target_compile_options(${module_name} PRIVATE ${BLAS_LINKER_FLAGS})
    target_link_libraries(${module_name} INTERFACE ${externalLib}::${externalLib})
else()
    set(${externalLib}_FOUND OFF)
endif (${externalLib}_FOUND)
set (${externalLib}_FOUND ${externalLib}_FOUND PARENT_SCOPE)