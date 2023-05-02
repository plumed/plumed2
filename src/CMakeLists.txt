cmake_minimum_required(VERSION 3.20.2)

project(Plumed2 LANGUAGES CXX)

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(alldirs ${CMAKE_CURRENT_SOURCE_DIR})# CACHE INTERNAL FORCE)
set(dirs ${alldirs})# CACHE INTERNAL FORCE)
#message("${alldirs}")
#message("${dirs}")
#these are the diretory to compile that do not have a module.type inside
set(coredirs main config wrapper)
set(mayBeExternal blas lapack)
list(REMOVE_ITEM dirs lib include maketools CMakeFiles ${coredirs} ${mayBeExternal})


foreach(dir ${coredirs})
    set(module_${dir} ON CACHE INTERNAL "always active core module ${dir}")
endforeach(dir ${coredirs})

foreach(dir ${mayBeExternal})
    option(module_${dir} "activate module ${dir}" ON)
endforeach(dir ${mayBeExternal})

option(all_modules "Activates all modules, if on ignores the values of module_name when compiling" OFF)

foreach(dir ${dirs})
    #this may the set up to be done in personal subdiretctories
    file(STRINGS ${dir}/module.type mt )
    if(mt STREQUAL "default-off" AND NOT all_modules)
        option(module_${dir} "activate module ${dir}" OFF)
    elseif(mt STREQUAL "default-on" OR (mt STREQUAL "default-off" AND all_modules))
        option(module_${dir} "activate module ${dir}" ON)
    elseif(mt STREQUAL "always")
         set(module_${dir} ON CACHE INTERNAL "always active module ${dir}")
    else() 
        #message( "${dir} ${mt}")
        message(FATAL_ERROR "module:\"${dir}\" type:\"${mt}\": module type must be \"always\" or \"deafult-on\" or \"default-off\"")
    endif()
    #this gets the dependency for the modules
    #file(STRINGS ${dir}/Makefile rgx REGEX "USE=.*")
    #STRING(REPLACE "USE=" "" rgx "${rgx}")
    #separate_arguments(rgx)
endforeach(dir ${dirs})
if(all_modules)
#forceactivates all modules
    foreach(dir ${dirs})
        set(module_${dir} ON)
    endforeach(dir ${dirs})
endif()
#  --enable-external-molfile-plugins
#                          enable search for external molfile_plugins, default:
#                          yes

option(useExternal_BLAS "enable search for external BLAS, default ON" ON)
option(useExternal_LAPACK "enable search for external LAPACK, default ON" ON)

foreach(externalLib BLAS LAPACK)
    string(TOLOWER ${externalLib} dir)
    if (useExternal_${externalLib})
        include(Find${externalLib})
    else()
        set(${externalLib}_FOUND OFF)
        set(module_${dir} ON CACHE BOOL "activate module ${dir}" FORCE)
    endif(useExternal_${externalLib})
    if (${externalLib}_FOUND)
        set(module_${dir} OFF CACHE BOOL "activate module ${dir}" FORCE)
        #target_compile_options(${PROJECT_NAME} PRIVATE ${BLAS_LINKER_FLAGS})
        #target_link_libraries(${PROJECT_NAME} ${BLAS_LIBRARIES})
endif ()
endforeach(externalLib BLAS LAPACK)


message("BLAS internal:${module_blas}")
message("LAPACK internal:${module_lapack}")
#options ## options are already done
##  --bindir=DIR            user executables [EPREFIX/bin]
##  --sbindir=DIR           system admin executables [EPREFIX/sbin]
# --libexecdir=DIR        program executables [EPREFIX/libexec]
##  --sysconfdir=DIR        read-only single-machine data [PREFIX/etc]
#  --sharedstatedir=DIR    modifiable architecture-independent data [PREFIX/com]
#  --localstatedir=DIR     modifiable single-machine data [PREFIX/var]
##  --libdir=DIR            object code libraries [EPREFIX/lib]
##  --includedir=DIR        C header files [PREFIX/include]
#  --oldincludedir=DIR     C header files for non-gcc [/usr/include]
#  --datarootdir=DIR       read-only arch.-independent data root [PREFIX/share]
#  --datadir=DIR           read-only architecture-independent data [DATAROOTDIR]
#  --infodir=DIR           info documentation [DATAROOTDIR/info]
#  --localedir=DIR         locale-dependent data [DATAROOTDIR/locale]
#  --mandir=DIR            man documentation [DATAROOTDIR/man]
#  --docdir=DIR            documentation root [DATAROOTDIR/doc/plumed]
#  --htmldir=DIR           html documentation [DOCDIR]
#  --dvidir=DIR            dvi documentation [DOCDIR]
#  --pdfdir=DIR            pdf documentation [DOCDIR]
#  --psdir=DIR             ps documentation [DOCDIR]
#
#Program names:
#  --program-prefix=PREFIX            prepend PREFIX to installed program names
#  --program-suffix=SUFFIX            append SUFFIX to installed program names
#  --program-transform-name=PROGRAM   run sed PROGRAM on installed program names
#
#Optional Features:
#  --disable-option-checking  ignore unrecognized --enable/--with options
#  --disable-FEATURE       do not include FEATURE (same as --enable-FEATURE=no)
#  --enable-FEATURE[=ARG]  include FEATURE [ARG=yes]
##  --enable-modules        all/none/reset or : separated list such as
##                          +crystallization:-bias default: reset
#  --enable-libsearch      enable search for libraries, default: yes
#  --enable-static-patch   enable allow statically linking plumed, default: yes
#  --enable-doc            enable documentation, default: yes
#  --enable-pdfdoc         enable pdf version of the manual, default: no
#  --enable-debug          enable debugging, default: no
#  --enable-gcov           enable gcov to estimate code coverage, default: no
#  --enable-cxx            11 or 14. To link libraries with headers that need
#                          this C++ level. Use --enable-cxx=none to remove
#                          -std=c++ flag
#  --enable-basic-warnings enable basic warnings, default: yes
#  --enable-fussy          enable fussy warnings, default: no
#  --enable-debug-glibcxx  enable enable boundary check, default: no
#  --enable-shared         enable shared libs, default: yes
#  --enable-dependency-tracking
#                          enable dependency tracking, default: yes
#  --enable-rpath          enable store rpath, default: no
#  --enable-absolute-soname
#                          enable store absolute soname (Linux only - this is
#                          the default behavior on OSX). Only enable for
#                          testing!, default: no
#  --enable-absolute-install-name
#                          enable store absolute relative (OSX only - disable
#                          to have a behavior similar to Linux). Only disable
#                          for testing!, default: yes
#  --enable-loader-path    enable use @loader_path to find
#                          libplumedKernel.dylib (OSX only), default: yes
#  --enable-bsymbolic      enable use -Bsymbolic flag in making shared
#                          libraries (Linux only), default: yes
#  --enable-ld-r           enable group object files, default: yes
#  --enable-ar-cr          enable use ar to build libplumedWrapper.a, default:
#                          yes
#  --enable-static-archive enable try to build libplumed.a for static linking,
#                          default: yes
#  --enable-asmjit         enable enable embedded asmjit, default: yes
#  --enable-mpi            enable search for mpi, default: yes
##  --enable-external-lapack
##                          enable search for external lapack, default: yes
##  --enable-external-blas  enable search for external blas, default: yes
#  --enable-molfile-plugins
#                          enable use molfile_plugins, default: yes
#  --enable-external-molfile-plugins
#                          enable search for external molfile_plugins, default:
#                          yes
#  --enable-zlib           enable search for zlib, default: yes
#  --enable-readdir-r      enable search for readdir_r (threadsafe), default:
#                          no
#  --enable-cregex         enable search for C regular expression, default: yes
#  --enable-dlopen         enable search for dlopen, default: yes
#  --enable-rtld_default   enable search for RTLD_DEFAULT macro, default: yes
#  --enable-chdir          enable search for chdir function, default: yes
#  --enable-subprocess     enable search for functions needed to manage a
#                          subprocess, default: yes
#  --enable-getcwd         enable search for getcwd function, default: yes
#  --enable-popen          enable search for popen, default: yes
#  --enable-execinfo       enable search for execinfo, default: yes
#  --enable-gsl            enable search for gsl, default: yes
#  --enable-boost_graph    enable search for boost graph, default: no
#  --enable-boost_serialization
#                          enable search for boost serialization, default: no
#  --enable-fftw           enable search for fftw, default: yes
#  --enable-python         enable search for python, default: yes
#  --enable-af_ocl         enable search for arrayfire_ocl, default: no
#  --enable-af_cuda        enable search for arrayfire_cuda, default: no
#  --enable-af_cpu         enable search for arrayfire_cpu, default: no
#  --enable-libtorch       enable search for libtorch, default: no
#  --disable-openmp        do not use OpenMP

#for debugging purpose
#foreach(dir ${dirs} ${coredirs} ${mayBeExternal})
#    message("${dir}: ${module_${dir}}")
#endforeach(dir ${dirs} ${coredirs} ${mayBeExternal})