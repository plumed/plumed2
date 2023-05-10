#this is not a module, this is set up to make easier the dependency interface
set(module_name main)
set(module_${module_name} ON CACHE INTERNAL "always active core module ${module_name}")

#dynamically linked plumed
#$(PLUMED_MAIN_SHARED): $(OBJ_MAIN) $(OBJ_WRAPPER) $(PLUMED_KERNEL)
add_executable(plumed_bin main.cpp)#OBJ_MAIN
target_link_libraries(plumed_bin PUBLIC
    sharedplumedKernel #PLUMED_KERNEL
    PlumedStaticWrapper #OBJ_WRAPPER
    )
set_target_properties(plumed_bin
    PROPERTIES RUNTIME_OUTPUT_NAME ${program_name}
    )
target_include_directories(plumed_bin PRIVATE ${PLUMED_SRC})
#statically linked plumed
#$(PLUMED_MAIN_STATIC): $(OBJ_MAIN) $(OBJ_KERNEL) $(OBJ_WRAPPER)
add_executable(plumed_bin_static main.cpp) #$(OBJ_MAIN)
target_link_libraries(plumed_bin_static PUBLIC
    PlumedStaticWrapper #$(OBJ_WRAPPER)
    archiveplumedKernel
    #${modulesForKernel} #$(OBJ_KERNEL) 
    #Config #extra options (are within OBJ_KERNEL in original Makefile)
)
set_target_properties(plumed_bin_static
    PROPERTIES RUNTIME_OUTPUT_NAME ${program_name}-static
    )
target_include_directories(plumed_bin_static PRIVATE ${PLUMED_SRC})
#runtime linked plumed
#$(PLUMED_MAIN_RUNTIME): $(OBJ_MAIN) $(OBJ_DYNAMIC_WRAPPER)
add_executable(plumed_bin_runtime main.cpp)
target_link_libraries(plumed_bin_runtime
    PUBLIC
    PlumedDynamicWrapper #OBJ_DYNAMIC_WRAPPER
    )
set_target_properties(plumed_bin_runtime
    PROPERTIES RUNTIME_OUTPUT_NAME ${program_name}-runtime
    )
target_include_directories(plumed_bin_runtime PRIVATE ${PLUMED_SRC})

#consider adding install manifest:
#install(CODE "string(REPLACE \";\" \"\\n\" MY_CMAKE_INSTALL_MANIFEST_CONTENT \"\$\{CMAKE_INSTALL_MANIFEST_FILES\}\")\n\
#              file(WRITE my_install_manifest.txt \"\$\{MY_CMAKE_INSTALL_MANIFEST_CONTENT\}\")")
#thanks https://stackoverflow.com/a/66804151
