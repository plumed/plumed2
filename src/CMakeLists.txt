cmake_minimum_required(VERSION 3.20.2)

project(Plumed2 LANGUAGES C CXX Fortran)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS None Debug Release RelWithDebInfo MinSizeRel)
#some settingsvalues
set(program_name "plumed" CACHE STRING "the name of the main executable")
set(SOEXT "so" CACHE STRING "The extension of dynamic libraries (so/dylib)")
set_property(CACHE SOEXT PROPERTY STRINGS "so" "dylib")

set (PLUMED_SRC ${CMAKE_CURRENT_SOURCE_DIR})
set (PLUMED_MAKETOOLS ${PLUMED_SRC}/maketools )
cmake_path (GET PLUMED_SRC PARENT_PATH PLUMED_MAIN_DIR)
set (PLUMED_SCRIPTS ${PLUMED_MAIN_DIR}/scripts)
set (PLUMED_PATCHES ${PLUMED_MAIN_DIR}/patches)

option(all_modules "Activates all modules, if on ignores the values of module_name when compiling" OFF)

set(modulesForKernel "")

################################################################################
#Setting up global options
################################################################################

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
#Macro definitions
################################################################################
MACRO(SUBDIRLIST result curdir)
  #This macro is from a lot of posts on stackoverflow
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

MACRO(ADDMODULETOKERNEL module_name)
    #use: ADDMODULETOKERNEL(module_name listOfSources)
    #Please write the source files explicitly
    if(${module_${module_name}} )
        set(OTHERARGS ${ARGV})
        list(REMOVE_ITEM OTHERARGS ${module_name})
        add_library(${module_name} OBJECT ${OTHERARGS} )
        target_include_directories(${module_name} PRIVATE ${PLUMED_SRC})
        list(APPEND modulesForKernel ${module_name})
        set(modulesForKernel ${modulesForKernel} PARENT_SCOPE)
    endif(${module_${module_name}})
ENDMACRO(ADDMODULETOKERNEL)

MACRO(ADDMODULENEEDS module_name)
    #use: ADDMODULENEEDS(module_name listOfModules)

    set(OTHERARGS ${ARGV})
    list(REMOVE_ITEM OTHERARGS ${module_name})
    set(moduleNeeds_${module_name} ${OTHERARGS} PARENT_SCOPE)
    
ENDMACRO(ADDMODULENEEDS)

MACRO(ADDMODULEDEPENDENCIES module_name)
    #use: ADDMODULEDEPENDENCIES(module_name listOfModules)
    #adds the module on which this depends on (for compile definitions and so on)
    #NB you may create circular dependencies
    #NB this is particularly experimental
    if(${module_${module_name}} )
        set(OTHERARGS ${ARGV})
        list(REMOVE_ITEM OTHERARGS ${module_name})
        foreach(lib ${OTHERARGS})
            #message("${module_name} is linked with ${lib}")
            target_link_libraries(${module_name} PUBLIC ${lib})
        endforeach(lib ${OTHERARGS})
    endif(${module_${module_name}})
ENDMACRO(ADDMODULEDEPENDENCIES)

function(print_target_property target_name property)
    get_target_property(tempvar ${target_name} ${property})
    message("${target_name} <${property}>: ${tempvar}")
    unset(tempvar)    
endfunction(print_target_property)

################################################################################
#Setting up Modules
################################################################################

SUBDIRLIST(alldirs ${CMAKE_CURRENT_SOURCE_DIR})# CACHE INTERNAL FORCE)
set(dirs ${alldirs})# CACHE INTERNAL FORCE)
#directories that do not contain modules
set(utildirs lib include maketools cmake)
#these are the diretory to compile that do not have a module.type inside
set(execDirs main config wrapper)
#these are directory with codes that may be external
set(mayBeExternal blas lapack molfile)
list(REMOVE_ITEM dirs ${utildirs} ${execDirs} ${mayBeExternal})

foreach(dir ${execDirs})
    add_subdirectory(${dir})
endforeach(dir ${execDirs})

message("Checking if it is possible to find external alternatives for modules:")
#needed for FindMOLFILE
set (CMAKE_MODULE_PATH ${PLUMED_SRC}/cmake)
foreach(dir ${mayBeExternal})
    add_subdirectory(${dir})
endforeach(dir ${mayBeExternal})

foreach(dir ${dirs})
    add_subdirectory(${dir})
endforeach(dir ${dirs})
if(all_modules)
    message("All module will be forced ON")
    foreach(dir ${dirs})
        set(module_${dir} ON)
    endforeach(dir ${dirs})
endif()

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
#foreach(dir ${dirs} ${execDirs} ${mayBeExternal})
#    message("${dir}: ${module_${dir}}")
#endforeach(dir ${dirs} ${execDirs} ${mayBeExternal})

#check module dependencies
foreach(dir ${dirs} ${mayBeExternal})
    if(${module_${dir}})
        set(notFound "")
        #message("${dir} ${moduleNeeds_${dir}}")
        foreach(otherModule ${moduleNeeds_${dir}})
            if (NOT module_${otherModule})
                list(FIND mayBeExternal ${otherModule} isExternal)
                if(${isExternal})
                    string(TOUPPER ${otherModule} capsName)
                    if(${capsName}_FOUND)
                        continue()
                    endif(${capsName}_FOUND)
                endif(${isExternal})
                list(APPEND notFound ${otherModule})
                set (modulesDependenciesProblems ON)
            endif(NOT module_${otherModule})
        endforeach(otherModule moduleNeeds_${dir})
        if (notFound)
            message(SEND_ERROR "Module \"${dir}\" needs the following modules active: ${notFound}")
        endif(notFound)
        unset(notFound)
    endif(${module_${dir}})
endforeach(dir ${dirs})
if (${modulesDependenciesProblems})
    message(FATAL_ERROR "Problems in module dependencies")
endif(${modulesDependenciesProblems})
################################################################################
#Setting up the main libraries
################################################################################

set(KernelTargets "")

foreach(activeModule ${modulesForKernel})
    list(APPEND KernelTargets $<TARGET_OBJECTS:${activeModule}>)
endforeach(activeModule ${modulesForKernel})

#these targets build the two libPlumedKernel.so
#building the shared libraries with the collected modules
add_library(libplumedKernel SHARED)
set_target_properties(libplumedKernel
    PROPERTIES 
        LIBRARY_OUTPUT_NAME ${program_name}Kernel
        SUFFIX .${SOEXT}
    )

target_link_libraries(libplumedKernel PUBLIC ${modulesForKernel})

add_library(libplumed SHARED $<TARGET_OBJECTS:PlumedWrapper>)
set_target_properties(libplumed
    PROPERTIES
        LIBRARY_OUTPUT_NAME ${program_name}
        SUFFIX .${SOEXT}
    )
target_link_libraries(libplumed
    PlumedWrapper
    ConfigInstall
)

add_library(libplumedWrapper STATIC $<TARGET_OBJECTS:PlumedWrapper>)
set_target_properties(libplumedWrapper
    PROPERTIES
        LIBRARY_OUTPUT_NAME ${program_name}Wrapper
        ARCHIVE_OUTPUT_NAME ${program_name}Wrapper
    )
target_link_libraries(libplumedWrapper
    PlumedWrapper
    ConfigInstall
)

#print_target_property(cltools COMPILE_DEFINITIONS)
#print_target_property(molfile COMPILE_DEFINITIONS)

#TODO move this in a better place:
IF(${OPT_DEBUGWARNINGS})
    if(MSVC)
      add_compile_options(/W4 /WX)
    else()
      add_compile_options(
          -Wall
          -Wextra
          -Wshadow
          -Wnon-virtual-dtor
          -Wold-style-cast
          -Wcast-align
          -Wunused
          -Woverloaded-virtual
          -Wpedantic
          -Wconversion
          -Wsign-conversion
          -Wnull-dereference
          -Wdouble-promotion
          -Wformat=2
          -Wduplicated-cond
          -Wduplicated-branches
          -Wlogical-op
          -Wuseless-cast
          #-Wlifetime
          )
    endif()
endif()
#Warnings
#-Wall
#-Wextra # reasonable and standard
#-Wshadow # warn the user if a variable declaration shadows one from a# parent context
#-Wnon-virtual-dtor # warn the user if a class with virtual functions has a# non-virtual destructor.
#-Wold-style-cast # warn for c-style casts
#-Wcast-align # warn for potential performance problem casts
#-Wunused # warn on anything being unused
#-Woverloaded-virtual # warn if you overload (not override) a virtual func
#-Wpedantic # warn if non-standard C++ is used
#-Wconversion # warn on type conversions that may lose data
#-Wsign-conversion # warn on sign conversions
#-Wnull-dereference # warn if a null dereference is detected
#-Wdouble-promotion # warn if float is implicit promoted to double
#-Wformat=2 # warn on security issues around functions that format output# (ie printf)
#-Wduplicated-cond # warn if if / else chain has duplicated conditions
#-Wduplicated-branches # warn if if / else branches have duplicated code
#-Wlogical-op # warn about logical operations being used where bitwise were# probably wanted
#-Wuseless-cast # warn if you perform a cast to the same type
#-Wlifetime # ///

# foreach(module_name blas lapack)
# print_target_property(${module_name} COMPILE_OPTIONS)
# print_target_property(${module_name} COMPILE_DEFINITIONS)
# print_target_property(${module_name} INTERFACE_COMPILE_DEFINITIONS)
# endforeach()
################################################################################
#One can also use -DCMAKE_EXPORT_COMPILE_COMMANDS=ON in that case, such that we 
#have all the commands in the file compile_commands.json without actually 
#building the sources.
#-DCMAKE_EXPORT_COMPILE_COMMANDS=ON
#thanks https://stackoverflow.com/questions/2670121/using-cmake-with-gnu-make-how-can-i-see-the-exact-commands#comment112625585_2673355
################################################################################
