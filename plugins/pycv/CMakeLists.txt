cmake_minimum_required(VERSION 3.15...3.27)
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

message(
  STATUS
    "Everithing should work fine if you are in the same environment in which you have compiled plumed"
)
# FinPlumed is here:
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
# Finding necessary packages
find_package(Python 3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG REQUIRED)
find_package(Plumed REQUIRED)

# Finding optionals things
if(Plumed_HAS_MPI)
  find_package(MPI REQUIRED)
endif()
if(MPI_CXX_FOUND)
  list(APPEND extraLibs MPI::MPI_CXX)
endif()

if(Plumed_HAS_OPENMP)
  find_package(OpenMP REQUIRED)
endif()
if(OpenMP_CXX_FOUND)
  list(APPEND extraLibs OpenMP::OpenMP_CXX)
endif()

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fno-gnu-unique USE_NO_GNU_UNIQUE)
if(USE_NO_GNU_UNIQUE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-gnu-unique")
endif()

# plumed_STATIC_LDFLAGS_OTHER:INTERNAL=-rdynamic;-Wl,-Bsymbolic;-fopenmp
# -rdynamic is automatically set by cmake, and also -fPIC

################################################################################
################################the pycv library################################
################################################################################

add_library(
  PythonCVInterface SHARED src/ActionWithPython.cpp src/PythonCVInterface.cpp
                           src/PythonFunction.cpp)
set(PLUMED_MODULE_IDENTITY "" CACHE FILEPATH "Generated module identity source")
set(PLUMED_RUNTIME_IDENTITY "" CACHE FILEPATH "Generated Python runtime identity source")
# public, so they trickle down to the python module
target_compile_definitions(PythonCVInterface PUBLIC ${Plumed_CFLAGS} __PLUMED_HAS_PYTHON_EMBED=1)
target_include_directories(PythonCVInterface PUBLIC src ${Plumed_INCLUDEDIR})
# ####################################################################
# uncommenting this brings problems since some symbols here are needed
# by the python module even if it should be the correct setting...
# https://gcc.gnu.org/wiki/Visibility could be a starting point
#######################################################################
# target_compile_options(PythonCVInterface PRIVATE -fvisibility=hidden)
if(APPLE AND PLUMED_MODULE_IDENTITY)
  # The runtime helper is loaded first via the module dependency manifest.
  # Python APIs then resolve to an existing host interpreter, if present.
  add_library(PythonCVRuntime SHARED src/PythonRuntime.cpp "${PLUMED_RUNTIME_IDENTITY}")
  target_include_directories(PythonCVRuntime PRIVATE ${Plumed_INCLUDEDIR} "${Plumed_INCLUDEDIR}/plumed")
  target_link_libraries(PythonCVRuntime PRIVATE pybind11::embed plumedKernel)
  set_target_properties(PythonCVRuntime PROPERTIES PREFIX "" INSTALL_RPATH "@loader_path/../../..")
  install(TARGETS PythonCVRuntime DESTINATION pycv)
  target_link_libraries(PythonCVInterface PRIVATE pybind11::headers PythonCVRuntime)
  target_include_directories(PythonCVInterface PRIVATE ${Python_INCLUDE_DIRS})
  target_link_options(PythonCVInterface PRIVATE "-Wl,-undefined,dynamic_lookup")
  target_compile_definitions(PythonCVInterface PRIVATE PLUMED_PYCV_RUNTIME_HELPER=1)
else()
  target_link_libraries(PythonCVInterface PRIVATE pybind11::embed)
endif()
target_link_libraries(PythonCVInterface PUBLIC plumedKernel ${extraLibs})
# This is a plugin library, with the native shared-library suffix.
set_target_properties(PythonCVInterface PROPERTIES PREFIX "")

# Optional registration for installations managed as PLUMED modules.
if(PLUMED_MODULE_IDENTITY)
  target_sources(PythonCVInterface PRIVATE "${PLUMED_MODULE_IDENTITY}")
  target_include_directories(PythonCVInterface PRIVATE "${Plumed_INCLUDEDIR}/plumed")
endif()

install(TARGETS PythonCVInterface DESTINATION pycv)

################################################################################
###########################The pvCV companion module############################
################################################################################

pybind11_add_module(plumedCommunications src/PlumedPythonEmbeddedModule.cpp)
target_link_libraries(plumedCommunications PRIVATE pybind11::headers)
target_link_libraries(plumedCommunications PUBLIC PythonCVInterface)

# The install directory is the output (wheel) directory
install(TARGETS plumedCommunications DESTINATION .)

# Keep the companion and its plugin relocatable together. The companion sits
# in site-packages and the plugin in its pycv/ subdirectory.
if(APPLE)
  set_target_properties(PythonCVInterface PROPERTIES INSTALL_RPATH "@loader_path;@loader_path/../../..")
  set_target_properties(plumedCommunications PROPERTIES INSTALL_RPATH "@loader_path/pycv;@loader_path/../..")
else()
  set_target_properties(PythonCVInterface PROPERTIES INSTALL_RPATH "$ORIGIN/../../..")
  set_target_properties(plumedCommunications PROPERTIES INSTALL_RPATH "$ORIGIN/pycv;$ORIGIN/../..")
endif()
