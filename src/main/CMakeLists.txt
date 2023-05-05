set(module_name main)
#this is not a module, this is set up to make easier the dependency interface
set(module_${module_name} ON CACHE INTERNAL "always active core module ${module_name}")
add_executable(plumed_bin main.cpp)
target_link_libraries(plumed_bin libplumedKernel libplumed)
set_target_properties(plumed_bin
    PROPERTIES RUNTIME_OUTPUT_NAME ${program_name}
    )
target_include_directories(plumed_bin PRIVATE ${PLUMED_SRC})
