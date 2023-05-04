#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "molfile")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" ON)
ADDMODULETOKERNEL(${module_name}
crdplugin.cpp
dcdplugin.cpp
gromacsplugin.cpp
pdbplugin.cpp
)
ADDMODULENEEDS(${module_name}

)

#The following lines should mimik the following line
#CPPFLAGS+=-DSTATIC_PLUGIN  -DVMDPLUGIN=molfile_$(patsubst %.o,%,$@)
target_compile_definitions(${module_name} PRIVATE STATIC_PLUGIN)
get_target_property(MY_PROJECT_SOURCES ${module_name} SOURCES)
foreach(sourcefile ${MY_PROJECT_SOURCES})
    string(REPLACE ".cpp" "" simpleName ${sourcefile})
    
    set_source_files_properties(${sourcefile} PROPERTIES COMPILE_DEFINITIONS VMDPLUGIN=molfile_${simpleName})
    
endforeach(sourcefile in ${MY_PROJECT_SOURCES})

#this should propagate the compile definition to the main
target_compile_definitions(${module_name} PUBLIC __PLUMED_HAS_MOLFILE_PLUGINS)