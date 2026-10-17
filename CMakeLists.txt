# PLUMED (master / 2.11-dev) - experimental MSVC/Windows native CMake build
# This is NOT the official build system (upstream uses autoconf+make, Unix only).
# It builds PLUMED as a STATIC kernel library (PlumedStatic wrapper, no dlopen)
# so it can be linked directly into GROMACS.
#
# Status: experimental, probing -> production.
# Modules: "always" + "default-on" sets, mirroring upstream module.type files.
cmake_minimum_required(VERSION 3.18)
project(PlumedMSVC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(PLUMED_SRC ${CMAKE_CURRENT_SOURCE_DIR})
set(PLUMED_GEN ${CMAKE_CURRENT_BINARY_DIR}/gen)
file(MAKE_DIRECTORY ${PLUMED_GEN})

# ---------------------------------------------------------------- modules
set(ALWAYS_MODULES asmjit blas core lapack lepton small_vector tools valtools)
set(DEFAULT_ON_MODULES bias cltools colvar function generic gridtools isdb
    liquid_crystal matrixtools molfile multicolvar refdist secondarystructure
    setup vatom xdrfile)
set(MODULES ${ALWAYS_MODULES} ${DEFAULT_ON_MODULES})
# Optional default-off modules (upstream: enabled by touching <mod>.on)
option(PLUMED_ENABLE_VES "Compile the VES module (variational enhanced sampling)" OFF)
if(PLUMED_ENABLE_VES)
    list(APPEND MODULES ves)
endif()

# ------------------------------------------------- build-time .inc generation
# Python replication of upstream maketools scripts (makecmd, makeModuleMap, xxd)
find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(GEN_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cmake-msvc/generate_incs.py)
set(PPM_INC ${PLUMED_GEN}/PlumedMainEnum.inc ${PLUMED_GEN}/PlumedMainMap.inc
            ${PLUMED_GEN}/GREXEnum.inc ${PLUMED_GEN}/GREXMap.inc
            ${PLUMED_GEN}/CLToolMainEnum.inc ${PLUMED_GEN}/CLToolMainMap.inc
            ${PLUMED_GEN}/ModuleMap.inc)
# Config.inc / version.h / Makefile.conf.xxd: upstream generates them with sed
# (src/config/Makefile) and update-version.sh; here the same substitutions are
# done by the generator script, into the build tree.
set(CONFIG_OUT ${PLUMED_GEN}/Config.inc ${PLUMED_GEN}/version.h
               ${PLUMED_GEN}/Makefile.conf.xxd)
add_custom_command(
    OUTPUT ${PPM_INC} ${CONFIG_OUT}
    COMMAND ${Python3_EXECUTABLE} ${GEN_SCRIPT}
            --src ${PLUMED_SRC} --outdir ${PLUMED_GEN}
    DEPENDS ${GEN_SCRIPT}
            ${PLUMED_SRC}/src/core/PlumedMain.cpp
            ${PLUMED_SRC}/src/core/GREX.cpp
            ${PLUMED_SRC}/src/core/CLToolMain.cpp
            ${PLUMED_SRC}/src/config/Config.inc.in
            ${PLUMED_SRC}/VERSION.txt
    COMMENT "Generating PLUMED enum/map/config files"
    VERBATIM)

# completion.xxd: shell completion bytes (generated once from completion.sh;
# content is inert data for the Completion cltool)
set(COMPLETION_XXD ${PLUMED_GEN}/completion.xxd)
add_custom_command(
    OUTPUT ${COMPLETION_XXD}
    COMMAND ${Python3_EXECUTABLE} ${GEN_SCRIPT}
            --src ${PLUMED_SRC} --outdir ${PLUMED_GEN} --completion-only
    DEPENDS ${GEN_SCRIPT} ${PLUMED_SRC}/src/cltools/completion.sh
    COMMENT "Generating completion.xxd"
    VERBATIM)

# Config.inc / version.h / Makefile.conf.xxd are produced by the command above
# (see CONFIG_OUT); no generated files are kept in the source tree.

# ---------------------------------------------------------------- library
set(PLUMED_SOURCES)
set(PLUMED_INCLUDE_DIRS ${PLUMED_SRC}/src)
foreach(m ${MODULES} config wrapper main)
    aux_source_directory(${PLUMED_SRC}/src/${m} _srcs)
    list(APPEND PLUMED_SOURCES ${_srcs})
    list(APPEND PLUMED_INCLUDE_DIRS ${PLUMED_SRC}/src/${m})
endforeach()
# static wrapper only (no dlopen); drop the dynamic .c wrapper
list(FILTER PLUMED_SOURCES EXCLUDE REGEX "wrapper/Plumed\\.c$")
# Upstream links EITHER Config.o (non-install build) OR ConfigInstall.o
# (install build, see src/lib/Makefile: Config$(PLUMED_INSTALL).o). We build
# the non-install kernel, so ConfigInstall.cpp must not be linked at all --
# it defines the same functions with install-time paths baked in.
list(FILTER PLUMED_SOURCES EXCLUDE REGEX "config/ConfigInstall\\.cpp$")
# main.cpp defines main() and belongs to the plumed executable only. Keep it
# OUT of the kernel library so that MSVC /WHOLEARCHIVE of the library is
# conflict-free for MD-engine hosts (GROMACS provides its own main()).
list(FILTER PLUMED_SOURCES EXCLUDE REGEX "main/main\\.cpp$")

add_library(plumedKernel STATIC ${PLUMED_SOURCES} ${PPM_INC} ${COMPLETION_XXD} ${CONFIG_OUT})
target_compile_features(plumedKernel PUBLIC cxx_std_17)

target_compile_definitions(plumedKernel PRIVATE
    NDEBUG
    # lepton/windowsIncludes.h defaults to dllimport on MSVC when no mode is
    # defined, which breaks static linking; select the static-lib mode:
    LEPTON_BUILDING_STATIC_LIBRARY
    # optional features off for this build:
    # (no __PLUMED_HAS_DLOPEN, __PLUMED_HAS_DLADDR, __PLUMED_HAS_MPI,
    #  __PLUMED_HAS_ZLIB, __PLUMED_HAS_MATHEVAL, __PLUMED_HAS_MOLFILE_PLUGINS)
)
# include dirs: src root (for "core/X.h" style), each module dir (for
# "tools/Tools.h" resolved from sibling modules), and gen dir for the .inc files
target_include_directories(plumedKernel PRIVATE
    ${PLUMED_SRC}/src
    ${PLUMED_GEN}
    ${PLUMED_SRC}/src/config
    ${PLUMED_SRC}/src/main
    ${PLUMED_SRC}/src/wrapper)
foreach(m ${MODULES})
    target_include_directories(plumedKernel PRIVATE ${PLUMED_SRC}/src/${m})
endforeach()

if(MSVC)
    target_compile_options(plumedKernel PRIVATE
        /permissive-      # enable and/or/not alternative tokens
        /EHsc
        /utf-8
        /Zc:__cplusplus   # MSVC reports real __cplusplus; Plumed.h relies on it
        /W1)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(plumedKernel PRIVATE /O2 /MD)
    endif()
endif()

# ---------------------------------------------------------------- main tool
add_executable(plumed ${PLUMED_SRC}/src/main/main.cpp)
# Registration objects in a static lib are dropped unless referenced:
# every RegisterMe global object is unreferenced by design, so pull in the
# whole library (upstream links all kernel objects directly, not an archive).
target_link_options(plumed PRIVATE
    "/WHOLEARCHIVE:$<TARGET_FILE:plumedKernel>")
target_link_libraries(plumed PRIVATE plumedKernel)
target_include_directories(plumed PRIVATE
    ${PLUMED_SRC}/src ${PLUMED_GEN} ${PLUMED_SRC}/src/config)
if(MSVC)
    target_compile_options(plumed PRIVATE /permissive- /EHsc /utf-8 /Zc:__cplusplus /W1)
endif()

list(LENGTH PLUMED_SOURCES PLUMED_NSRC)
message(STATUS "PLUMED MSVC build: ${PLUMED_NSRC} sources")
