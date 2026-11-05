cmake_minimum_required(VERSION 3.16)
project(plumed_thread_mpi LANGUAGES CXX)

find_package(Threads REQUIRED)
# Imported upstream implementation sources stay in sync with reimport.py.
file(GLOB TMPI_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_library(plumed_thread_mpi STATIC ${TMPI_SOURCES})
target_include_directories(plumed_thread_mpi PUBLIC . PRIVATE .)
target_compile_features(plumed_thread_mpi PUBLIC cxx_std_17)
target_compile_definitions(plumed_thread_mpi PRIVATE
  THREAD_PTHREADS HAVE_UNISTD_H HAVE_SYS_TIME_H HAVE_PTHREAD_H
  SIZEOF_LONG_LONG_INT __PLUMED_HAS_THREAD_MPI)
target_link_libraries(plumed_thread_mpi PUBLIC Threads::Threads)

set(TEST_INPUTS
  launch.cpp
  mpi-bindings.cpp
  host-runner.cpp
  internal-runner.cpp
  internal-launcher.cpp
  isolation-main.cpp
  mpiwrap-runtime.cpp)
foreach(name IN LISTS TEST_INPUTS)
  configure_file(${name}.in ${name} COPYONLY)
endforeach()

add_executable(thread_mpi_launch_test ${CMAKE_CURRENT_BINARY_DIR}/launch.cpp)
target_link_libraries(thread_mpi_launch_test PRIVATE plumed_thread_mpi)

add_executable(thread_mpi_bindings_test
  ${CMAKE_CURRENT_BINARY_DIR}/mpi-bindings.cpp)
target_include_directories(thread_mpi_bindings_test PRIVATE .)
target_link_libraries(thread_mpi_bindings_test PRIVATE plumed_thread_mpi)

add_executable(thread_mpi_internal_launcher_test
  ${CMAKE_CURRENT_BINARY_DIR}/internal-launcher.cpp
  ../mpiwrap/mpiwrap.cpp
  ../mpiwrap/mpiwrap-launcher.cpp
  ../tools/Exception.cpp)
target_include_directories(thread_mpi_internal_launcher_test PRIVATE
  .. ../mpiwrap .)
target_compile_definitions(thread_mpi_internal_launcher_test PRIVATE
  SIZEOF_LONG_LONG_INT __PLUMED_HAS_THREAD_MPI __PLUMED_MPIWRAP)
target_link_libraries(thread_mpi_internal_launcher_test PRIVATE plumed_thread_mpi)

add_executable(mpiwrap_runtime_test
  ${CMAKE_CURRENT_BINARY_DIR}/mpiwrap-runtime.cpp
  ../mpiwrap/mpiwrap.cpp
  ../mpiwrap/mpiwrap-launcher.cpp
  ../tools/Exception.cpp)
target_include_directories(mpiwrap_runtime_test PRIVATE
  .. ../mpiwrap .)
target_compile_definitions(mpiwrap_runtime_test PRIVATE
  SIZEOF_LONG_LONG_INT __PLUMED_HAS_THREAD_MPI __PLUMED_MPIWRAP)
target_link_libraries(mpiwrap_runtime_test PRIVATE plumed_thread_mpi)

enable_testing()
foreach(ranks IN ITEMS 1 2 4)
  add_test(NAME launch_${ranks} COMMAND thread_mpi_launch_test ${ranks})
  set_tests_properties(launch_${ranks} PROPERTIES TIMEOUT 20)
  add_test(NAME bindings_${ranks} COMMAND thread_mpi_bindings_test ${ranks})
  set_tests_properties(bindings_${ranks} PROPERTIES TIMEOUT 20)
  add_test(NAME internal_launcher_${ranks}
    COMMAND thread_mpi_internal_launcher_test ${ranks})
  set_tests_properties(internal_launcher_${ranks} PROPERTIES TIMEOUT 20)
endforeach()
add_test(NAME internal_launcher_callback_exception
  COMMAND ${CMAKE_COMMAND}
  -DLAUNCHER=$<TARGET_FILE:thread_mpi_internal_launcher_test>
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check-launcher-fatal.cmake)
set_tests_properties(internal_launcher_callback_exception PROPERTIES
  TIMEOUT 20)
add_test(NAME internal_unlaunched_serial
  COMMAND thread_mpi_internal_launcher_test serial)
set_tests_properties(internal_unlaunched_serial PROPERTIES TIMEOUT 20)
foreach(ranks IN ITEMS 1 2 4)
  add_test(NAME mpiwrap_runtime_${ranks} COMMAND mpiwrap_runtime_test ${ranks})
  set_tests_properties(mpiwrap_runtime_${ranks} PROPERTIES TIMEOUT 20)
endforeach()

# Compile the host and serial kernel separately: only the host links tMPI.
# This exercises the generated Plumed.h descriptor across a DSO boundary,
# rather than installing the embedded provider in the same translation unit.
if(UNIX)
  configure_file(injected-smoke-host.cpp.in injected-smoke-host.cpp COPYONLY)
  configure_file(injected-smoke-kernel.cpp.in injected-smoke-kernel.cpp COPYONLY)
  add_library(mpiwrap_injected_smoke_kernel SHARED
    ${CMAKE_CURRENT_BINARY_DIR}/injected-smoke-kernel.cpp ../mpiwrap/mpiwrap.cpp)
  target_include_directories(mpiwrap_injected_smoke_kernel PRIVATE .. ../mpiwrap)
  target_compile_features(mpiwrap_injected_smoke_kernel PRIVATE cxx_std_17)
  target_compile_definitions(mpiwrap_injected_smoke_kernel PRIVATE __PLUMED_MPIWRAP)
  target_link_libraries(mpiwrap_injected_smoke_kernel PRIVATE Threads::Threads)
  add_executable(mpiwrap_injected_smoke_host
    ${CMAKE_CURRENT_BINARY_DIR}/injected-smoke-host.cpp)
  target_include_directories(mpiwrap_injected_smoke_host PRIVATE .. .)
  target_compile_definitions(mpiwrap_injected_smoke_host PRIVATE SIZEOF_LONG_LONG_INT)
  target_link_libraries(mpiwrap_injected_smoke_host PRIVATE
    plumed_thread_mpi ${CMAKE_DL_LIBS})
  add_test(NAME injected_thread_mpi
    COMMAND mpiwrap_injected_smoke_host $<TARGET_FILE:mpiwrap_injected_smoke_kernel>)
  set_tests_properties(injected_thread_mpi PROPERTIES TIMEOUT 20
    ENVIRONMENT "PLUMED_MPI_WRAPPER=;PLUMED_MPI_WRAPPER_FORCE=")
endif()
