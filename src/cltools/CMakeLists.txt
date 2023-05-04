#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "cltools")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" ON)
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
#generated file
completion.xxd
)
ADDMODULEDEPENDENCIES(${module_name}
core config tools molfile xdrfile
)

#get_target_property(t ${module_name} INCLUDE_DIRECTORIES)
#message("id: ${t}")
target_include_directories(${module_name} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#get_target_property(t ${module_name} INCLUDE_DIRECTORIES)
#message("id: ${t}")
#this file once compiled reproduces
#completion.xxd: completion.sh ../../scripts/*.sh ../../patches/*.sh
#	{ ../maketools/make-scripts-options.sh ; cat completion.sh ; }| ../maketools/xxd  > completion.xxd
configure_file(completion.xxd.compiler.sh.in completion.xxd.compiler.sh @ONLY)
#this is used to test completion.xxd.compiler.sh
#execute_process(
#        COMMAND bash completion.xxd.compiler.sh
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#        #OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/completion.xxd
#)

file(GLOB script_deps ${PLUMED_SCRIPTS}/*.sh)
file(GLOB patches_deps ${PLUMED_PATCHES}/*.sh)

add_custom_command(OUTPUT completion.xxd
    DEPENDS completion.sh ${script_deps} ${patches_deps}
    COMMAND bash completion.xxd.compiler.sh
    COMMENT "Generationc completion.xxd"
    #WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

#message("scripts ${script_deps} ${PLUMED_SCRIPTS}")
#message("patches ${patches_deps} ${PLUMED_PATCHES}")
