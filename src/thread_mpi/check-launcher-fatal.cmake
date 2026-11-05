execute_process(
  COMMAND "${LAUNCHER}" 2 throw
  RESULT_VARIABLE launch_status
  OUTPUT_VARIABLE launch_stdout
  ERROR_VARIABLE launch_stderr
  TIMEOUT 10)

if(launch_status STREQUAL "0")
  message(FATAL_ERROR "rank callback exception unexpectedly returned success")
endif()
if(NOT launch_stderr MATCHES "PLUMED thread-MPI rank 1 failed: injected callback failure")
  message(FATAL_ERROR "rank callback failure diagnostic missing: ${launch_stderr}")
endif()
