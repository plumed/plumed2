#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "blas")
string(TOUPPER ${module_name} externalLib)
#Note that the macros here require this directory added as a subdir of plumed/src
set(module_${module_name} ON CACHE INTERNAL "always active module ${module_name}")
option(useExternal_${externalLib} "enable search for external ${externalLib}, default ON" ON)
ADDMODULETOKERNEL(${module_name}
blas.cpp
)
#ADDMODULENEEDS(${module_name}
#
#)

#include(FortranCInterface)
if (useExternal_${externalLib})
    find_package(${externalLib})
endif (useExternal_${externalLib})

if (${externalLib}_FOUND)
    target_compile_definitions(${module_name}
        PUBLIC
        __PLUMED_HAS_EXTERNAL_BLAS=1)
    include(CheckSourceCompiles)
    #From the manual:The check is only performed once, with the result cached in
    #the variable named by <resultVar>. Every subsequent CMake run will re-use
    #this cached value rather than performing the check again, even if the
    #<code> changes. In order to force the check to be re-evaluated, the
    #variable named by <resultVar> must be manually removed from the cache.
    set(CMAKE_REQUIRED_LIBRARIES "BLAS::BLAS")
    check_source_compiles(C #using C to not use extern "C"
        "int main(void) {
            double a,b,c,s;
            srotg_(&a,&b,&c,&s);
        }"
    UnderscodeBlas)
    check_source_compiles(C  #using C to not use extern "C"
        "int main(void) {
            double a,b,c,s;
            srotg(&a,&b,&c,&s);
        }"
    noUnderscodeBlas
    )
    unset(CMAKE_REQUIRED_LIBRARIES)

    if(noUnderscodeBlas AND NOT UnderscodeBlas)
        target_compile_definitions(${module_name}
            PUBLIC
            F77_NO_UNDERSCORE)
    elseif(NOT noUnderscodeBlas AND NOT UnderscodeBlas)
        message(WARNING "both srtog_ and srtog are not linkable")
    elseif(noUnderscodeBlas AND UnderscodeBlas)
        message(WARNING "both srtog_ and srtog are linkable")
    endif(noUnderscodeBlas AND NOT UnderscodeBlas)

    target_link_libraries(${module_name}
        PUBLIC
        ${externalLib}::${externalLib})
else()
    set(${externalLib}_FOUND OFF)
endif (${externalLib}_FOUND)
set (${externalLib}_FOUND ${externalLib}_FOUND PARENT_SCOPE)

# print_target_property(${module_name} COMPILE_OPTIONS)
# print_target_property(${module_name} COMPILE_DEFINITIONS)
# print_target_property(${module_name} INTERFACE_COMPILE_DEFINITIONS)
