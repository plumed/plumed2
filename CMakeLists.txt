cmake_minimum_required(VERSION 3.20.2)



project(Plumed2 LANGUAGES CXX)

set(PLUMED_VERSION_SHORT "2.10")
set(PLUMED_VERSION_LONG "2.10.0-dev")
#set(PLUMED_VERSION_GIT "23e754fd4-dirty")
execute_process(
        COMMAND git describe --long --dirty --always
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        OUTPUT_VARIABLE PLUMED_VERSION_GIT
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(PLUMED_VERSION_MAJOR 2)
set(PLUMED_VERSION_MINOR 10)
set(PLUMED_VERSION_PATCH 0)
#as now CMake is only experimental, it can be useful if you use vscode or other
#IDE that have an interface with it, but we are not still confident that it is ready
#for installation purposes
set(CMAKE_SKIP_INSTALL_RULES YES)
option(risk_accepted "CMake is experimental use at your own risk" OFF)
if(NOT risk_accepted)
   message(FATAL_ERROR "You must accept the risk to use the experimental cmake compilation: -Drisk_accepted=ON")
endif(NOT risk_accepted)


add_subdirectory(src)
