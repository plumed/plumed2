#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "blas")
#Note that the macros here require this directory added as a subdir of plumed/src
set(module_${module_name} ON CACHE INTERNAL "always active module ${module_name}")
ADDMODULETOKERNEL(${module_name}
blas.cpp
)
ADDMODULENEEDS(${module_name}

)

#include(FortranCInterface)

string(TOUPPER ${module_name} externalLib)
if (useExternal_${externalLib})
    find_package(${externalLib})
endif (useExternal_${externalLib})

if (${externalLib}_FOUND)
    target_compile_definitions(${module_name} PUBLIC __PLUMED_HAS_EXTERNAL_BLAS=1)
    #activate/deactivate
#    (underscore) search_for=dsyevr_ ;;
#    (nounderscore) search_for=dsyevr ;;
    target_compile_definitions(${module_name} PUBLIC F77_NO_UNDERSCORE)
    target_compile_options(${module_name} PRIVATE ${BLAS_LINKER_FLAGS})
    target_link_libraries(${module_name} PUBLIC ${BLAS_LIBRARIES}) #${externalLib}::${externalLib})
else()
    set(${externalLib}_FOUND OFF)
endif (${externalLib}_FOUND)
set (${externalLib}_FOUND ${externalLib}_FOUND PARENT_SCOPE)