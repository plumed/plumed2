
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

add_library(Config OBJECT Config.cpp)
add_library(ConfigInstall OBJECT ConfigInstall.cpp)