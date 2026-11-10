# PLUMED (master / 2.11-dev) - experimental MSVC/Windows native CMake build
# This is NOT the official build system (upstream uses autoconf+make, Unix only).
# It builds PLUMED as a STATIC kernel library (PlumedStatic wrapper, no dlopen)
# so it can be linked directly into GROMACS.
#
# Status: experimental, probing -> production.
# Modules: "always" + "default-on" sets, mirroring upstream module.type files.
cmake_minimum_required(VERSION 3.18)
project(PlumedMSVC LANGUAGES CXX)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(PLUMED_SRC ${CMAKE_CURRENT_SOURCE_DIR})
set(PLUMED_GEN ${CMAKE_CURRENT_BINARY_DIR}/gen)
file(MAKE_DIRECTORY ${PLUMED_GEN})

# ---------------------------------------------------------------- modules
# The "always" and "default-on" modules, as listed by the src/*/module.type
# files (the same set that ./configure enables by default)
set(MODULES)
file(GLOB _module_type_files ${PLUMED_SRC}/src/*/module.type)
list(SORT _module_type_files)
foreach(_mt ${_module_type_files})
    file(STRINGS ${_mt} _type LIMIT_COUNT 1)
    if(_type STREQUAL "always" OR _type STREQUAL "default-on")
        get_filename_component(_dir ${_mt} DIRECTORY)
        get_filename_component(_mod ${_dir} NAME)
        list(APPEND MODULES ${_mod})
    endif()
endforeach()
# The Windows CMake build uses the imported Win32 thread backend by default.
# Keep the wrapper frontend enabled when this option is off so the fully
# serial build still provides the communicator compatibility layer.
option(PLUMED_ENABLE_THREAD_MPI "Compile the embedded thread-MPI backend" ON)
if(NOT PLUMED_ENABLE_THREAD_MPI)
    list(REMOVE_ITEM MODULES thread_mpi)
endif()
# Optional default-off modules (upstream: enabled by touching <mod>.on)
option(PLUMED_ENABLE_VES "Compile the VES module (variational enhanced sampling)" OFF)
if(PLUMED_ENABLE_VES)
    list(APPEND MODULES ves)
endif()

# ------------------------------------------------- build-time .inc generation
# CMake reimplementation of the upstream maketools scripts (makecmd,
# makeModuleMap, xxd, update-version.sh and the Config.inc rules from
# src/config/Makefile), so that no external interpreter (python/awk/sed) is
# needed at build time.
set(GEN_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cmake-msvc/generate.cmake)
set(GEN_MODULE ${CMAKE_CURRENT_SOURCE_DIR}/cmake-msvc/PlumedGenerate.cmake)
set(PPM_INC ${PLUMED_GEN}/PlumedMainEnum.inc ${PLUMED_GEN}/PlumedMainMap.inc
            ${PLUMED_GEN}/GREXEnum.inc ${PLUMED_GEN}/GREXMap.inc
            ${PLUMED_GEN}/CLToolMainEnum.inc ${PLUMED_GEN}/CLToolMainMap.inc
            ${PLUMED_GEN}/ModuleMap.inc
            ${PLUMED_GEN}/ActionDescriptions.inc)
set(CONFIG_OUT ${PLUMED_GEN}/Config.inc ${PLUMED_GEN}/version.h
               ${PLUMED_GEN}/Makefile.conf.xxd)
set(COMPLETION_XXD ${PLUMED_GEN}/completion.xxd)
# The module map is built from every module source; make the generation
# command depend on them so that it is rerun when actions are added/removed.
file(GLOB _generator_inputs CONFIGURE_DEPENDS ${PLUMED_SRC}/src/*/*.cpp)
add_custom_command(
    OUTPUT ${PPM_INC} ${CONFIG_OUT} ${COMPLETION_XXD}
    COMMAND ${CMAKE_COMMAND}
            "-DSRC=${PLUMED_SRC}" "-DOUT=${PLUMED_GEN}"
            -P ${GEN_SCRIPT}
    DEPENDS ${GEN_SCRIPT} ${GEN_MODULE}
            ${PLUMED_SRC}/src/config/Config.inc.in
            ${PLUMED_SRC}/VERSION.txt
            ${PLUMED_SRC}/src/cltools/completion.sh
            ${PLUMED_SRC}/src/maketools/makeActionDescriptions.cmake
            ${_generator_inputs}
    COMMAND ${CMAKE_COMMAND} "-DSRC=${PLUMED_SRC}" "-DOUT=${PLUMED_GEN}"
            -P ${PLUMED_SRC}/src/maketools/makeActionDescriptions.cmake
    COMMENT "Generating PLUMED enum/map/config files (CMake)"
    VERBATIM)

# ---------------------------------------------------------------- library
set(PLUMED_SOURCES)
set(PLUMED_INCLUDE_DIRS ${PLUMED_SRC}/src)
foreach(m ${MODULES} config wrapper main)
    # aux_source_directory appends to the variable
    set(_srcs)
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

# Shared kernel (plumedKernel.dll on Windows), loadable with PLUMED_KERNEL and
# plumed_create_dlopen (e.g. by the python interface). Only its C entry points are
# exported (cmake-msvc/plumedKernel.def); the plumed executable is linked to it.
option(PLUMED_SHARED_KERNEL "Build plumedKernel as a shared library and link plumed to it" OFF)
if(PLUMED_SHARED_KERNEL)
    set(_kernel_type SHARED)
else()
    set(_kernel_type STATIC)
endif()
add_library(plumedKernel ${_kernel_type} ${PLUMED_SOURCES} ${PPM_INC} ${COMPLETION_XXD} ${CONFIG_OUT})
if(PLUMED_SHARED_KERNEL AND MSVC)
    target_sources(plumedKernel PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cmake-msvc/plumedKernel.def)
endif()
target_compile_features(plumedKernel PUBLIC cxx_std_17)

target_compile_definitions(plumedKernel PRIVATE
    NDEBUG
    __PLUMED_MPIWRAP=1
    # lepton/windowsIncludes.h defaults to dllimport on MSVC when no mode is
    # defined, which breaks static linking; select the static-lib mode:
    LEPTON_BUILDING_STATIC_LIBRARY
    # optional features off for this build:
    # (no __PLUMED_HAS_DLOPEN, __PLUMED_HAS_DLADDR, __PLUMED_HAS_MPI,
    #  __PLUMED_HAS_ZLIB, __PLUMED_HAS_MATHEVAL)
)
# Embedded asmjit, used by lepton to compile expressions (CUSTOM, MATHEVAL,
# switching functions...) to machine code. As in ./configure, only on x86-64
# and on 64 bit ARM; PLUMED_USE_ASMJIT=no disables it at run time.
option(PLUMED_ENABLE_ASMJIT "Compile expressions with the embedded asmjit" ON)
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _processor)
if(PLUMED_ENABLE_ASMJIT AND _processor MATCHES "^(x86_64|amd64|aarch64|arm64)$"
   AND CMAKE_SIZEOF_VOID_P EQUAL 8)
    target_compile_definitions(plumedKernel PRIVATE __PLUMED_HAS_ASMJIT)
    if(UNIX AND NOT APPLE)
        # clock_gettime, shm_open and shm_unlink are in librt on older Linux systems
        find_library(PLUMED_RT_LIBRARY rt)
        if(PLUMED_RT_LIBRARY)
            target_link_libraries(plumedKernel PRIVATE ${PLUMED_RT_LIBRARY})
        endif()
    endif()
    set(PLUMED_HAS_ASMJIT ON)
else()
    set(PLUMED_HAS_ASMJIT OFF)
endif()
message(STATUS "PLUMED asmjit: ${PLUMED_HAS_ASMJIT} (${CMAKE_SYSTEM_PROCESSOR})")
if(PLUMED_ENABLE_THREAD_MPI)
    target_compile_definitions(plumedKernel PRIVATE
        __PLUMED_HAS_THREAD_MPI=1
        SIZEOF_LONG_LONG_INT
    )
    if(WIN32)
        target_compile_definitions(plumedKernel PRIVATE THREAD_WINDOWS)
    elseif(UNIX)
        target_compile_definitions(plumedKernel PRIVATE
            THREAD_PTHREADS HAVE_UNISTD_H HAVE_SYS_TIME_H HAVE_PTHREAD_H
        )
        find_package(Threads REQUIRED)
        target_link_libraries(plumedKernel PUBLIC Threads::Threads)
    else()
        message(FATAL_ERROR "PLUMED_ENABLE_THREAD_MPI requires a supported thread backend")
    endif()
endif()
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
        /W1
        /wd4068)          # unknown pragma (#pragma acc); same as -Wno-unknown-pragmas in configure
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(plumedKernel PRIVATE /O2 /MD)
    endif()
endif()
# debug information in an optimized build, to get stack traces of crashes
# (CI: plumed.exe run under cdb). /Z7 keeps it in the objects of the static
# kernel library, the linker collects it in plumed.pdb.
option(PLUMED_MSVC_DEBUG_INFO "MSVC: debug information also in Release builds" OFF)
if(MSVC AND PLUMED_MSVC_DEBUG_INFO)
    target_compile_options(plumedKernel PRIVATE /Z7)
    if(PLUMED_SHARED_KERNEL)
        target_link_options(plumedKernel PRIVATE /DEBUG /OPT:REF /OPT:ICF)
    endif()
endif()

# ------------------------------------------------------------ config.txt
# Read by "plumed config" (scripts/config.sh) from PLUMED_ROOT/src/config:
# regtests use it to check for features (plumed_needs) and modules
# (plumed_modules). Same format as src/maketools/update-config-txt.sh,
# restricted to what this build can enable.
set(_has cregex)
if(PLUMED_HAS_ASMJIT)
    list(APPEND _has asmjit)
endif()
if(PLUMED_ENABLE_THREAD_MPI)
    list(APPEND _has mpi thread_mpi)
endif()
set(_config_txt "# generated by the CMake build (CMakeLists.txt)\n\n")
string(APPEND _config_txt "python_bin \n\nmpiexec \n\n")
foreach(_h ${_has})
    string(APPEND _config_txt "has ${_h} on\n")
endforeach()
string(APPEND _config_txt "\n")
file(GLOB _module_types ${PLUMED_SRC}/src/*/module.type)
list(SORT _module_types)
foreach(_mt ${_module_types})
    get_filename_component(_dir ${_mt} DIRECTORY)
    get_filename_component(_mod ${_dir} NAME)
    file(STRINGS ${_mt} _type LIMIT_COUNT 1)
    if(_mod IN_LIST MODULES)
        string(APPEND _config_txt "module ${_mod} on (${_type})\n")
    else()
        string(APPEND _config_txt "module ${_mod} off (${_type})\n")
    endif()
endforeach()
file(CONFIGURE OUTPUT ${PLUMED_SRC}/src/config/config.txt CONTENT "${_config_txt}" @ONLY)

# ---------------------------------------------------------------- main tool
if(PLUMED_SHARED_KERNEL)
    # main.cpp uses the C interface only: give it its own wrapper, which gets
    # the kernel functions from the shared kernel (see cmake-msvc/PlumedLinked.cpp)
    add_executable(plumed ${PLUMED_SRC}/src/main/main.cpp
                          ${CMAKE_CURRENT_SOURCE_DIR}/cmake-msvc/PlumedLinked.cpp)
    if(APPLE)
        set_target_properties(plumed PROPERTIES INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
    elseif(UNIX)
        set_target_properties(plumed PROPERTIES INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
    endif()
else()
    add_executable(plumed ${PLUMED_SRC}/src/main/main.cpp)
endif()
if(MSVC)
    # UTF-8 as the code page of plumed.exe (file names on the command line)
    target_sources(plumed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cmake-msvc/plumed.manifest)
endif()
# Registration objects in a static lib are dropped unless referenced:
# every RegisterMe global object is unreferenced by design, so pull in the
# whole library (upstream links all kernel objects directly, not an archive).
# A shared kernel contains all its objects anyway.
if(PLUMED_SHARED_KERNEL)
    if(MSVC)
        target_link_options(plumed PRIVATE "/STACK:8388608")
        if(PLUMED_MSVC_DEBUG_INFO)
            target_compile_options(plumed PRIVATE /Z7)
            target_link_options(plumed PRIVATE /DEBUG /OPT:REF /OPT:ICF)
        endif()
    endif()
elseif(MSVC)
    target_link_options(plumed PRIVATE "/WHOLEARCHIVE:$<TARGET_FILE:plumedKernel>")
    # 8 MB main-thread stack, as on Linux (the MSVC default is 1 MB): nested
    # cmd() calls (DEBUG logRequestedAtoms, the PLUMED action) overflowed it
    target_link_options(plumed PRIVATE "/STACK:8388608")
    if(PLUMED_MSVC_DEBUG_INFO)
        target_compile_options(plumed PRIVATE /Z7)
        target_link_options(plumed PRIVATE /DEBUG /OPT:REF /OPT:ICF)
    endif()
elseif(APPLE)
    target_link_options(plumed PRIVATE "LINKER:-force_load,$<TARGET_FILE:plumedKernel>")
else()
    target_link_options(plumed PRIVATE
        "LINKER:--whole-archive,$<TARGET_FILE:plumedKernel>,--no-whole-archive")
endif()
target_link_libraries(plumed PRIVATE plumedKernel)
target_include_directories(plumed PRIVATE
    ${PLUMED_SRC}/src ${PLUMED_GEN} ${PLUMED_SRC}/src/config)
if(MSVC)
    target_compile_options(plumed PRIVATE /permissive- /EHsc /utf-8 /Zc:__cplusplus /W1 /wd4068)
endif()

set(_compilable ${PLUMED_SOURCES})
list(FILTER _compilable INCLUDE REGEX "\.(cpp|c|cxx)$")
list(LENGTH _compilable PLUMED_NSRC)
message(STATUS "PLUMED MSVC build: ${PLUMED_NSRC} sources")

# ---------------------------------------------------------------- install
# Mirror the essential parts of `make install` (src/lib/Makefile) for the
# static-kernel build. See cmake-msvc/PlumedInstall.cmake for the details.
include(GNUInstallDirs)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake-msvc/PlumedInstall.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake-msvc/PlumedGenerate.cmake)
plumed_read_version("${PLUMED_SRC}")

set(PLUMED_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
set(PLUMED_INSTALL_PROGRAM_DIR "${CMAKE_INSTALL_LIBDIR}/plumed")

install(TARGETS plumed RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if(PLUMED_SHARED_KERNEL)
    # Windows: plumedKernel.dll next to plumed.exe, import library plumedKernel.lib in lib
    install(TARGETS plumedKernel
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${PLUMED_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${PLUMED_INSTALL_LIBDIR})
else()
    install(FILES $<TARGET_FILE:plumedKernel> DESTINATION ${PLUMED_INSTALL_LIBDIR})
endif()

# headers, with the same include rewrite as maketools/install.headers
install(CODE "
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake-msvc/PlumedInstall.cmake)
  message(STATUS \"Installing headers\")
  plumed_install_headers(\"${PLUMED_SRC}\" \"\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/plumed\")
")

# link files for MD engines (make and cmake), pkg-config file, patches/scripts
set(_libplumeddir "\${CMAKE_INSTALL_PREFIX}/${PLUMED_INSTALL_PROGRAM_DIR}")
set(_srcroot "${CMAKE_CURRENT_SOURCE_DIR}")
set(_libdir "\${CMAKE_INSTALL_PREFIX}/${PLUMED_INSTALL_LIBDIR}")
install(CODE "
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake-msvc/PlumedInstall.cmake)
  file(MAKE_DIRECTORY \"${_libplumeddir}/src/lib\" \"\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/pkgconfig\")
  message(STATUS \"Writing Plumed.inc / Plumed.cmake\")
  plumed_write_link_files(\"${_libdir}/plumedKernel.lib\" \"${_libdir}\" plumed \"${_libplumeddir}/src/lib\" \"${_srcroot}\")
  message(STATUS \"Writing pkg-config file\")
  plumed_write_pkgconfig(plumed \"\${CMAKE_INSTALL_PREFIX}\" \"\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}\" \"${_libdir}\" \"${PLUMED_VERSION_LONG}\" \"\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/pkgconfig/plumed.pc\")
")
install(DIRECTORY ${PLUMED_SRC}/patches ${PLUMED_SRC}/scripts DESTINATION ${PLUMED_INSTALL_PROGRAM_DIR})
install(FILES ${PLUMED_SRC}/fortran/plumed.f90 ${PLUMED_SRC}/fortran/plumed_f08.f90 DESTINATION ${PLUMED_INSTALL_PROGRAM_DIR}/fortran)
