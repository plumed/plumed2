cmake_minimum_required(VERSION 3.0)
project(astyle CXX)

# Release Build - release by default (except for Borland)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Build Options - executable by default, libraries on request
option(BUILD_JAVA_LIBS   "Build java library"   OFF)
option(BUILD_SHARED_LIBS "Build shared library" OFF)
option(BUILD_STATIC_LIBS "Build static library" OFF)

# Linux Soname Version
set(MAJORVER 3)
set(MINORVER 1)
set(PATCHVER 0)
set(SOLIBVER ${MAJORVER}.${MINORVER}.${PATCHVER})

# AStyle Source
list(APPEND SRCS
    src/ASBeautifier.cpp
    src/ASEnhancer.cpp
    src/ASFormatter.cpp
    src/ASLocalizer.cpp
    src/ASResource.cpp
    src/astyle_main.cpp)

# AStyle Documentation
list(APPEND DOCS
    doc/astyle.html
    doc/install.html
    doc/news.html
    doc/notes.html
    doc/styles.css)

# If a java library is requested, shared libraries should be enabled
# and the Java Development Kit 'include' directories added
if(BUILD_JAVA_LIBS)
    set(BUILD_SHARED_LIBS ON)
    if(WIN32)
        set(java_home $ENV{JAVA_HOME})
        if(NOT java_home)
            message(FATAL_ERROR "Environment variable JAVA_HOME not defined")
        endif()
        if(NOT EXISTS ${java_home}/include)
            message(FATAL_ERROR "Java Development Kit not installed at ${java_home}")
        endif()
    else()
        if(NOT EXISTS /usr/lib/jvm/default-java/include)
            message(FATAL_ERROR "Java Development Kit not installed")
        endif()
    endif()
endif()

# Define the output type and install directories
# To uninstall 'xargs rm < install_manifest.txt'
if(BUILD_SHARED_LIBS OR BUILD_STATIC_LIBS)
    add_library(astyle ${SRCS})
    if(NOT WIN32)
        install(TARGETS astyle DESTINATION /usr/lib)
    endif()
else()
    add_executable(astyle ${SRCS})
    if(WIN32)
        set(pf86 "PROGRAMFILES(x86)")
        set(prog_files $ENV{${pf86}})
        if(NOT prog_files)
            set(prog_files $ENV{PROGRAMFILES})
        endif()
        install(TARGETS astyle DESTINATION ${prog_files}/AStyle)
        install(FILES ${DOCS} DESTINATION ${prog_files}/AStyle/doc/)
    else()
        install(TARGETS astyle DESTINATION /usr/bin)
        install(FILES ${DOCS} DESTINATION /usr/share/doc/astyle)
    endif()
endif()

# Set build-specific compile options
if(BUILD_SHARED_LIBS OR BUILD_STATIC_LIBS)
    if(BUILD_JAVA_LIBS)
        target_compile_options(astyle PRIVATE -DASTYLE_JNI)
        if(WIN32)
            target_include_directories(astyle PRIVATE $ENV{JAVA_HOME}/include)
            target_include_directories(astyle PRIVATE $ENV{JAVA_HOME}/include/win32)
        else()
            target_include_directories(astyle PRIVATE /usr/lib/jvm/default-java/include)
            target_include_directories(astyle PRIVATE /usr/lib/jvm/default-java/include/linux)
        endif()
    else()
        target_compile_options(astyle PRIVATE -DASTYLE_LIB)
    endif()
    # Windows DLL exports removed
    set_property(TARGET astyle PROPERTY DEFINE_SYMBOL "")
    # Linux solib version added
    set_property(TARGET astyle PROPERTY VERSION ${SOLIBVER})
    set_property(TARGET astyle PROPERTY SOVERSION ${MAJORVER})
endif()

# Set file names if different than 'astyle'
if(BUILD_JAVA_LIBS)
    if(WIN32)
        set_property(TARGET astyle PROPERTY OUTPUT_NAME AStyle31j)
        set_property(TARGET astyle PROPERTY PREFIX "")
    else()
        set_property(TARGET astyle PROPERTY OUTPUT_NAME astylej)
    endif()
elseif(BUILD_SHARED_LIBS)
    if(WIN32)
        set_property(TARGET astyle PROPERTY OUTPUT_NAME AStyle31)
        set_property(TARGET astyle PROPERTY PREFIX "")
    endif()
elseif(BUILD_STATIC_LIBS)
    if(WIN32)
        set_property(TARGET astyle PROPERTY OUTPUT_NAME AStyleLib)
        set_property(TARGET astyle PROPERTY PREFIX "")
    endif()
else()
    if(WIN32)
        set_property(TARGET astyle PROPERTY OUTPUT_NAME AStyle)
    endif()
endif()

macro(set_linker_options strip_option)
    # Remove -rdynamic and add 'strip' to linker flags
    if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
        set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
        set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${strip_option}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${strip_option}")
    endif()
    # Shared library options
    if(BUILD_SHARED_LIBS)
        if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
            set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-intel")
        elseif(MINGW)
            # minGW dlls don't work
            # tdm-gcc dlls work with everything except python
            set(CMAKE_SHARED_LINKER_FLAGS
                "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--add-stdcall-alias -Wl,--dll")
        elseif(BORLAND)
            # use a development environment to compile Borland dlls
        endif()
    endif()
endmacro()

# Set default compile options for supported compilers
if(APPLE)
    target_compile_options(
        astyle PRIVATE -W -Wall -fno-rtti -fno-exceptions -std=c++11 -stdlib=libc++)
    set_linker_options("")
elseif(NOT WIN32)   # Linux
    target_compile_options(astyle PRIVATE -Wall -fno-rtti -fno-exceptions -std=c++11)
    set_linker_options("-s")
elseif(BORLAND)     # Release must be explicitely requested for Borland
    target_compile_options(astyle PRIVATE -w -x-)   # Cannot use no-rtti (-RT-)
    set_linker_options("")
elseif(MINGW)
    target_compile_options(astyle PRIVATE -Wall -Wextra -fno-rtti -fno-exceptions -std=c++0x)
    set_linker_options("-s")
endif()

# Display build information
if(BUILD_JAVA_LIBS)
    message("CMAKE_BUILD_TYPE is Java ${CMAKE_BUILD_TYPE} ${SOLIBVER}")
elseif(BUILD_SHARED_LIBS)
    message("CMAKE_BUILD_TYPE is Shared ${CMAKE_BUILD_TYPE} ${SOLIBVER}")
elseif(BUILD_STATIC_LIBS)
    message("CMAKE_BUILD_TYPE is Static ${CMAKE_BUILD_TYPE}")
else()
    message("CMAKE_BUILD_TYPE is Executable ${CMAKE_BUILD_TYPE}")
endif()
set(CMAKE_VERBOSE_MAKEFILE ON)
