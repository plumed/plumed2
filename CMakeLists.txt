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

add_subdirectory(src)
