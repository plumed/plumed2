add_executable(plumed_bin main.cpp)
target_link_libraries(plumed_bin libplumedKernel libplumed)
set_target_properties(plumed_bin
    PROPERTIES RUNTIME_OUTPUT_NAME ${program_name}
    )