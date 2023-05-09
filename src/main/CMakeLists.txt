#this is not a module, this is set up to make easier the dependency interface
set(module_name main)
set(module_${module_name} ON CACHE INTERNAL "always active core module ${module_name}")

#dynamically linked plumed
add_executable(plumed_bin main.cpp)
target_link_libraries(plumed_bin PUBLIC libplumedKernel libplumed)
set_target_properties(plumed_bin
    PROPERTIES RUNTIME_OUTPUT_NAME ${program_name}
    )
target_include_directories(plumed_bin PRIVATE ${PLUMED_SRC})
#statically linked plumed
add_executable(plumed_bin_static main.cpp)
target_link_libraries(plumed_bin_static PUBLIC PlumedStaticWrapper)
target_link_libraries(plumed_bin_static PUBLIC ${modulesForKernel})
target_link_libraries(plumed_bin_static PUBLIC ConfigInstall)
set_target_properties(plumed_bin_static
    PROPERTIES RUNTIME_OUTPUT_NAME ${program_name}-static
    )
target_include_directories(plumed_bin_static PRIVATE ${PLUMED_SRC})
#runtime linked plumed
# add_executable(plumed_bin_runtime main.cpp)
# target_link_libraries(plumed_bin_runtime PUBLIC libplumedWrapper libplumedStatic)
# set_target_properties(plumed_bin_runtime
    # PROPERTIES RUNTIME_OUTPUT_NAME ${program_name}-runtime
    # )
# target_include_directories(plumed_bin_runtime PRIVATE ${PLUMED_SRC})
