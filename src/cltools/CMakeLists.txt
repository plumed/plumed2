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

#add_custom_command(OUTPUT completion.xxd
execute_process(
        COMMAND echo "\{ ${PLUMED_MAKETOOLS}/make-scripts-options.sh ; cat completion.sh; \} \| ${PLUMED_MAKETOOLS}/xxd  > ${CMAKE_CURRENT_BINARY_DIR}/completion.xxd"
        COMMAND_ECHO STDOUT
)

message (FATAL_ERROR "cltools has a non standard Makefile (more than 4 lines) you need to modify the CMakeLists.txt!")

#completion.xxd: completion.sh ../../scripts/*.sh ../../patches/*.sh
#q	{ ../maketools/make-scripts-options.sh ; cat completion.sh ; }| ../maketools/xxd  > completion.xxd