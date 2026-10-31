# CMake-native reimplementation of the PLUMED build-time generators.
#
# Mirrors maketools/makecmd, makeModuleMap, xxd, update-version.sh and the
# Config.inc / Makefile.conf rules from src/config/Makefile, so that no
# external interpreter (python/awk/sed) is needed at build time.
#
# Note: CMake regular expressions do not support POSIX character classes
# ([[:space:]]); use explicit character classes such as [ \t\r\n] instead.

# --- makecmd -------------------------------------------------------------
# Extract the "case cmd_xxx:" entries from a source file and write
# <base>Enum.inc / <base>Map.inc (same content as maketools/makecmd).
function(plumed_makecmd cpp out_dir base)
  file(STRINGS "${cpp}" _lines REGEX "^[ \t]*case[ \t]+cmd_[A-Za-z0-9_]+:")
  set(_enum "")
  set(_map "")
  foreach(_l IN LISTS _lines)
    string(REGEX MATCH "cmd_[A-Za-z0-9_]+" _name "${_l}")
    if(_name)
      string(REGEX REPLACE "^cmd_" "" _short "${_name}")
      list(APPEND _enum "${_name}")
      list(APPEND _map "  {\"${_short}\",${_name}}")
    endif()
  endforeach()
  string(REPLACE ";" ",\n" _enum_txt "${_enum}")
  string(REPLACE ";" ",\n" _map_txt "${_map}")
  file(WRITE "${out_dir}/${base}Enum.inc" "${_enum_txt}\n")
  file(WRITE "${out_dir}/${base}Map.inc" "${_map_txt}\n")
  list(LENGTH _enum _n)
  message(STATUS "  ${base}: ${_n} commands")
endfunction()

# --- makeModuleMap -------------------------------------------------------
# Map every PLUMED_REGISTER_ACTION(name) to the module directory that
# contains it (same content as maketools/makeModuleMap).
function(plumed_make_module_map src_root out_file)
  set(_rows "")
  file(GLOB _mt_files "${src_root}/src/*/module.type")
  foreach(_mt IN LISTS _mt_files)
    get_filename_component(_moddir "${_mt}" DIRECTORY)
    get_filename_component(_mod "${_moddir}" NAME)
    file(GLOB _cpps "${_moddir}/*.cpp")
    foreach(_cpp IN LISTS _cpps)
      file(READ "${_cpp}" _text)
      string(REGEX MATCHALL "PLUMED_REGISTER_ACTION[ \t\r\n]*[(][ \t\r\n]*[0-9A-Za-z_:]+[ \t\r\n]*,[ \t\r\n]*\"([0-9A-Za-z_]+)\"" _ms "${_text}")
      foreach(_m IN LISTS _ms)
        string(REGEX REPLACE "^.*\"([0-9A-Za-z_]+)\"$" "\\1" _action "${_m}")
        list(APPEND _rows "  {\"${_action}\", \"${_mod}\"}")
      endforeach()
    endforeach()
  endforeach()
  string(REPLACE ";" ",\n" _txt "${_rows}")
  file(WRITE "${out_file}" "${_txt}\n")
  list(LENGTH _rows _n)
  message(STATUS "  ModuleMap: ${_n} actions")
endfunction()

# --- xxd -----------------------------------------------------------------
# Byte dump in the format produced by maketools/xxd (12 bytes per line).
function(plumed_xxd_file in_file out_file)
  file(READ "${in_file}" _hex HEX)
  string(LENGTH "${_hex}" _len)
  math(EXPR _n "${_len}/2")
  set(_out "")
  if(_n GREATER 0)
    math(EXPR _last "${_n}-1")
    foreach(_i RANGE 0 ${_last})
      math(EXPR _off "${_i}*2")
      string(SUBSTRING "${_hex}" ${_off} 2 _b)
      if(_i GREATER 0)
        math(EXPR _mod "${_i}%12")
        if(_mod EQUAL 0)
          string(APPEND _out ",\n ")
        else()
          string(APPEND _out ",")
        endif()
      endif()
      string(APPEND _out " 0x${_b}")
    endforeach()
  endif()
  string(APPEND _out "\n")
  file(WRITE "${out_file}" "${_out}")
endfunction()

# --- version -------------------------------------------------------------
# Same values as maketools/update-version.sh: VERSION.txt + git describe.
function(plumed_read_version src_root)
  file(STRINGS "${src_root}/VERSION.txt" _lines)
  set(_ver "")
  foreach(_l IN LISTS _lines)
    if(NOT _l MATCHES "^#")
      string(STRIP "${_l}" _l)
      if(NOT _l STREQUAL "")
        set(_ver "${_l}")
        break()
      endif()
    endif()
  endforeach()
  if(_ver MATCHES "^([0-9]+)[.]([0-9]+)([.]([0-9]+))?")
    set(_major "${CMAKE_MATCH_1}")
    set(_minor "${CMAKE_MATCH_2}")
    if(CMAKE_MATCH_4)
      set(_patch "${CMAKE_MATCH_4}")
    else()
      set(_patch "0")
    endif()
    set(_short "${_major}.${_minor}")
  else()
    set(_short "Unknown")
    set(_major "0")
    set(_minor "0")
    set(_patch "0")
  endif()
  execute_process(COMMAND git describe --long --dirty --always
                  WORKING_DIRECTORY "${src_root}"
                  OUTPUT_VARIABLE _git OUTPUT_STRIP_TRAILING_WHITESPACE
                  RESULT_VARIABLE _rc ERROR_QUIET)
  if(NOT _rc EQUAL 0 OR _git STREQUAL "")
    execute_process(COMMAND git rev-parse --short=12 HEAD
                    WORKING_DIRECTORY "${src_root}"
                    OUTPUT_VARIABLE _git OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE _rc2 ERROR_QUIET)
    if(NOT _rc2 EQUAL 0 OR _git STREQUAL "")
      set(_git "Unknown")
    endif()
  endif()
  set(PLUMED_VERSION_SHORT "${_short}" PARENT_SCOPE)
  set(PLUMED_VERSION_LONG "${_ver}" PARENT_SCOPE)
  set(PLUMED_VERSION_GIT "${_git}" PARENT_SCOPE)
  set(PLUMED_VERSION_MAJOR "${_major}" PARENT_SCOPE)
  set(PLUMED_VERSION_MINOR "${_minor}" PARENT_SCOPE)
  set(PLUMED_VERSION_PATCH "${_patch}" PARENT_SCOPE)
endfunction()

# --- Config.inc / version.h / Makefile.conf.xxd --------------------------
function(plumed_generate_config src_root out_dir)
  plumed_read_version("${src_root}")

  # Config.inc: the substitutions that src/config/Makefile does with sed for
  # the non-installed build. configure_file() is the CMake equivalent.
  get_filename_component(_root "${src_root}" ABSOLUTE)
  file(TO_CMAKE_PATH "${_root}" _root)
  set(PLUMED_ROOT "${_root}")
  set(SOEXT "dll")
  set(ISINSTALLED "false")
  set(htmldir "xxxxNAxxxx")
  set(includedir "xxxxNAxxxx")
  set(program_name "plumed")
  configure_file("${src_root}/src/config/Config.inc.in" "${out_dir}/Config.inc" @ONLY)

  # version.h (update-version.sh)
  configure_file("${CMAKE_CURRENT_LIST_DIR}/version.h.in" "${out_dir}/version.h" @ONLY)

  # Makefile.conf, embedded through Makefile.conf.xxd (src/config/Makefile)
  configure_file("${CMAKE_CURRENT_LIST_DIR}/Makefile.conf.in" "${out_dir}/Makefile.conf" @ONLY)
  plumed_xxd_file("${out_dir}/Makefile.conf" "${out_dir}/Makefile.conf.xxd")

  message(STATUS "  Config.inc / version.h / Makefile.conf.xxd (version ${PLUMED_VERSION_LONG}, git ${PLUMED_VERSION_GIT})")
endfunction()

# --- completion.xxd ------------------------------------------------------
# Same content as { make-scripts-options.sh ; cat completion.sh ; } | xxd.
# The scripts are run only if bash is available; otherwise a minimal set of
# global options is embedded (the completion cltool only prints this text).
function(plumed_generate_completion src_root out_dir)
  find_program(_bash bash)
  set(_content "")
  if(_bash)
    file(GLOB _scripts "${src_root}/scripts/*.sh")
    foreach(_s IN LISTS _scripts)
      get_filename_component(_name "${_s}" NAME_WE)
      string(REPLACE "-" "_" _name "${_name}")
      execute_process(COMMAND "${_bash}" -c "PLUMED_ROOT='${src_root}' bash '${_s}' --options"
                      WORKING_DIRECTORY "${src_root}"
                      OUTPUT_VARIABLE _opts ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
      if(NOT _opts STREQUAL "")
        string(APPEND _content "local cmd_keys_${_name}=\"${_opts}\"\n")
      endif()
    endforeach()
  endif()
  if(_content STREQUAL "")
    set(_content "local cmd_keys_plumed=\"--help --options --version\"\n")
  endif()
  file(READ "${src_root}/src/cltools/completion.sh" _comp)
  string(APPEND _content "${_comp}")
  set(_tmp "${out_dir}/completion.txt")
  file(WRITE "${_tmp}" "${_content}")
  plumed_xxd_file("${_tmp}" "${out_dir}/completion.xxd")
  file(REMOVE "${_tmp}")
  message(STATUS "  completion.xxd (bash: ${_bash})")
endfunction()
