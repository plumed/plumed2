# Driver for the PLUMED build-time generators (run with: cmake -P).
# Usage: cmake -DSRC=<plumed-source-root> -DOUT=<generation-dir> -P generate.cmake
if(NOT DEFINED SRC OR NOT DEFINED OUT)
  message(FATAL_ERROR "Usage: cmake -DSRC=<plumed-source-root> -DOUT=<generation-dir> -P generate.cmake")
endif()
include("${CMAKE_CURRENT_LIST_DIR}/PlumedGenerate.cmake")
file(MAKE_DIRECTORY "${OUT}")

plumed_makecmd("${SRC}/src/core/PlumedMain.cpp"   "${OUT}" "PlumedMain")
plumed_makecmd("${SRC}/src/core/GREX.cpp"         "${OUT}" "GREX")
plumed_makecmd("${SRC}/src/core/CLToolMain.cpp"   "${OUT}" "CLToolMain")
plumed_make_module_map("${SRC}" "${OUT}/ModuleMap.inc")
plumed_generate_config("${SRC}" "${OUT}")
plumed_generate_completion("${SRC}" "${OUT}")
