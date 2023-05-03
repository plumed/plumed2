#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "cltools")
#Note that the macros here require this directory added as a subdir of plumed/src

option(module_${module_name} "activate module cltools" ON)
ADDMODULETOKERNEL(${module_name}
CLTool.cpp
Completion.cpp
Driver.cpp
DriverDouble.cpp
DriverFloat.cpp
GenExample.cpp
GenJson.cpp
GenTemplate.cpp
Info.cpp
kT.cpp
Manual.cpp
PdbRenumber.cpp
pesmd.cpp
SimpleMD.cpp
SumHills.cpp
completion.xxd
)
ADDMODULEDEPENDENCIES(${module_name}
core config tools molfile xdrfile
)

#TODO: this is a workaround, need to be done as add_custom_command
#add_custom_command(OUTPUT completion.xxd
execute_process(
        COMMAND ${PLUMED_MAKETOOLS}/make-scripts-options.sh
        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/make_scripts_options
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
execute_process(
        COMMAND cat ${CMAKE_CURRENT_BINARY_DIR}/make_scripts_options ${CMAKE_CURRENT_SOURCE_DIR}/completion.sh    
        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/completions
)

execute_process(
        COMMAND ${PLUMED_MAKETOOLS}/xxd
        INPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/completions
        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/completion.xxd
)
#message("${cat_make_scripts_options} ${cat_completion}")
message (FATAL_ERROR "cltools has a non standard Makefile (more than 4 lines) you need to modify the CMakeLists.txt!")

#completion.xxd: completion.sh ../../scripts/*.sh ../../patches/*.sh
#q	{ ../maketools/make-scripts-options.sh ; cat completion.sh ; }| ../maketools/xxd  > completion.xxd