if(NOT Plumed_FOUND)
  find_package(PkgConfig)
  if(Plumed_FIND_QUIETLY)
    function(message)
      # THIS completely shuts down messages
    endfunction()
    pkg_check_modules(PLUMED QUIET plumedInternals)
  else()
    pkg_check_modules(PLUMED plumedInternals)
  endif()

  if(Plumed_FOUND)
    if("-D__PLUMED_HAS_MPI=1" IN_LIST Plumed_CFLAGS)
      set(Plumed_HAS_MPI
          1
          CACHE INTERNAL "plumed has MPI")
    endif()
    if("-fopenmp" IN_LIST Plumed_STATIC_LDFLAGS_OTHER)
      set(Plumed_HAS_OPENMP
          1
          CACHE INTERNAL "plumed has OpenMP")
    endif()
  else()
    message(STATUS "plumed not found via pkgconfig, trying executable")

    execute_process(
      COMMAND plumed --no-mpi info --include-dir
      RESULT_VARIABLE PLUMED_EXECUTABLE
      OUTPUT_QUIET ERROR_QUIET)
    if(PLUMED_EXECUTABLE EQUAL 0)
      set(Plumed_FOUND
          1
          CACHE INTERNAL "plumed found")

      message(STATUS "Configuring plumed from executable")

      execute_process(
        COMMAND plumed --no-mpi info --configuration
        OUTPUT_VARIABLE Plumed_CONFIG
        OUTPUT_STRIP_TRAILING_WHITESPACE)

      set(Plumed_CPP_FLAGS "")
      set(Plumed_EXTRA_INCLUDE_FLAGS "")
      set(tPlumed_CPP_FLAGS "")

      string(REPLACE "\n" ";" ProcessFile_LINES "${Plumed_CONFIG}")
      foreach(_line ${ProcessFile_LINES})
        if(${_line} MATCHES "CPPFLAGS=.*")
          set(tPlumed_CPP_FLAGS ${_line})
          string(REGEX REPLACE "CPPFLAGS= *" "" tPlumed_CPP_FLAGS
                               ${tPlumed_CPP_FLAGS})
          string(REPLACE "\\" "" tPlumed_CPP_FLAGS ${tPlumed_CPP_FLAGS})
          string(REPLACE " -" ";-" tPlumed_CPP_FLAGS ${tPlumed_CPP_FLAGS})

          # message(STATUS "Found PLUMED CPP_FLAGS: \'${tPlumed_CPP_FLAGS}\'")
          # message(STATUS "Found PLUMED CPP_FLAGS:")
          foreach(_flag ${tPlumed_CPP_FLAGS})
            # message(STATUS "   \"${_flag}\"")
            if(${_flag} MATCHES "^-D")
              string(REPLACE "-D" "" t ${_flag})
              list(APPEND Plumed_CPP_FLAGS ${t})
            endif()
            if(${_flag} MATCHES "^-I")
              string(REPLACE "-I" "" t ${_flag})
              list(APPEND Plumed_EXTRA_INCLUDE_FLAGS ${t})
            endif()
          endforeach()
        endif()
        if(${_line} MATCHES ".*-fopenmp.*")
          set(Plumed_HAS_MPI
              1
              CACHE INTERNAL "plumed has MPI")
        endif()
      endforeach()

      set(Plumed_CFLAGS
          ${Plumed_CPP_FLAGS}
          CACHE INTERNAL "plumed Definitions flags")

      execute_process(
        COMMAND plumed --no-mpi info --include-dir
        OUTPUT_VARIABLE Plumed_INCLUDEDIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)

      list(APPEND Plumed_INCLUDEDIR ${Plumed_EXTRA_INCLUDE_FLAGS})

      set(Plumed_INCLUDEDIR
          ${Plumed_INCLUDEDIR}
          CACHE INTERNAL "plumed include dir")

      execute_process(COMMAND plumed --no-mpi config -q has mpi
                      RESULT_VARIABLE Plumed_WITH_MPI)
      if(Plumed_WITH_MPI EQUAL 0)
        set(Plumed_HAS_MPI
            1
            CACHE INTERNAL "plumed has MPI")
      endif()

    else()
      if(Plumed_FIND_REQUIRED)
        message(FATAL_ERROR "plumed not found")
      endif()
    endif()
  endif()
endif()
