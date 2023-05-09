set(module_name config)
#this is not a module, this is set up to make easier the dependency interface
set(module_${module_name} ON CACHE INTERNAL "always active core module ${module_name}")
#SOEXT is set in ./src
set(ISINSTALLED "true")
#must go to /lib/plumed
set(PLUMED_ROOT ${CMAKE_INSTALL_PREFIX}/${program_name})
set (htmldir ${htmldir})
set (includedir ${includedir})
set (program_name ${program_name})
set (libdir ${libdir})
configure_file(Config.inc.in ConfigInstall.inc)

set(ISINSTALLED "false")
set (PLUMED_ROOT ${CMAKE_SOURCE_DIR})
set (htmldir "xxxxNAxxxx")
set (includedir "xxxxNAxxxx")
set (program_name "xxxxNAxxxx")
configure_file(Config.inc.in Config.inc)

#TODO: verify if the version number are not set
configure_file(version.h.in version.h)
################################################################################
#The config interface library holds a lot of options that must be inherited by #
#everithing, more specialized options will be found in the respective modules  #
# licke blas and lapack
################################################################################
add_library(config INTERFACE)
################################################################################
#                              options                                         #
################################################################################
option(enable_dlopen "enable search for dlopen" ON)
if(enable_dlopen)
    if(CMAKE_DL_LIBS)
        target_link_libraries(config INTERFACE ${CMAKE_DL_LIBS})
        target_compile_definitions(config INTERFACE __PLUMED_HAS_DLOPEN)
    endif(CMAKE_DL_LIBS)
endif(enable_dlopen)
#  --disable-openmp        do not use OpenMP
option(enable_openmp "enable search for openmp" ON)
if(enable_openmp)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(config INTERFACE OpenMP::OpenMP_CXX)
    endif()
endif(enable_openmp)

#must be set up a target interface config wotj all the definitios
add_library(Config OBJECT Config.cpp Makefile.conf.xxd)
target_link_libraries(Config PUBLIC config)
target_include_directories(Config PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

add_library(ConfigInstall OBJECT ConfigInstall.cpp Makefile.conf.xxd)
target_link_libraries(ConfigInstall PUBLIC config)
target_include_directories(ConfigInstall PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

configure_file(Makefile.conf.xxd.compiler.sh.in Makefile.conf.xxd.compiler.sh @ONLY)
configure_file(${PLUMED_MAIN_DIR}/Makefile.conf.in Makefile.cmakeMocked.conf @ONLY)
add_custom_command(OUTPUT Makefile.conf.xxd
    DEPENDS ${PLUMED_MAIN_DIR}/Makefile.conf.in
    COMMAND bash Makefile.conf.xxd.compiler.sh
    COMMENT "Generationc completion.xxd"
    #WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

#print_target_property(ConfigInstall INCLUDE_DIRECTORIES)
