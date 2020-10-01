# DO NOT EDIT ./configure FILE
# One should edit configure.ac and the run autoconf in this
# directory. Also the resulting "configure" is stored on the git,
# so as to allow people not to install a recent autoconf on their
# system.
# If you modify configure.ac, remember to run
# autoconf and then commit both files to the git repository.

# we require a recent version
# notice that autoconf is not necessary on user's machine, but only
# if one has to update configure.ac
AC_PREREQ([2.68])
AC_INIT([PLUMED], [2])

##################################################################
# In order to add m4 macros, put then in a m4 directory and include them here.
# E.g.:
# m4_include([./m4/m4_ax_openmp.m4])

##################################################################
# Here we define a few useful macros

# PLUMED_CONFIG_ENABLE(variablename,optionname,doc,default)
# notice that variablename and optionname are likely identical,
# they just need to be different with optionname contains a "-"
# (not allowed in shell variable names)
AC_DEFUN([PLUMED_CONFIG_ENABLE], [
m4_bpatsubst([$1],-,_)=
AC_ARG_ENABLE([$1],
  AS_HELP_STRING([--enable-$1], [enable $2, default: $3]),
  [case "${enableval}" in
             (yes) m4_bpatsubst([$1],-,_)=true ;;
             (no)  m4_bpatsubst([$1],-,_)=false ;;
             (*)   AC_MSG_ERROR([wrong argument to --enable-$1]) ;;
  esac],
  [case "$3" in
             (yes) m4_bpatsubst([$1],-,_)=true ;;
             (no)  m4_bpatsubst([$1],-,_)=false ;;
  esac]
)
])

# PLUMED_CHECK_CXXFLAG(flag)
# use it to check if a flag is available on this compiler
AC_DEFUN([PLUMED_CHECK_CXXFLAG], [
  save_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$CXXFLAGS $1"
  AC_MSG_CHECKING([whether $CXX accepts $1])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [
      AC_LINK_IFELSE([AC_LANG_PROGRAM([])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([not linking]); CXXFLAGS="$save_CXXFLAGS"])
    ],
    [AC_MSG_RESULT([no]); CXXFLAGS="$save_CXXFLAGS"]
  )
])

# PLUMED_CHECK_CFLAG(flag)
# use it to check if a flag is available on this compiler
AC_DEFUN([PLUMED_CHECK_CFLAG], [
  AC_LANG_PUSH(C)
  save_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS $1"
  AC_MSG_CHECKING([whether $CC accepts $1])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [
      AC_LINK_IFELSE([AC_LANG_PROGRAM([])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([not linking]); CCFLAGS="$save_CCFLAGS"])
    ],
    [AC_MSG_RESULT([no]); CCFLAGS="$save_CCFLAGS"]
  )
  AC_LANG_POP
])

# PLUMED_CHECK_LDFLAGS(flag)
# use it to check if a flag is available on this compiler
AC_DEFUN([PLUMED_CHECK_LDFLAGS], [
  save_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $1"
  AC_MSG_CHECKING([whether LDFLAGS can contain $1])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no]); LDFLAGS="$save_LDFLAGS"])
])

# PLUMED_SEARCH_LIBS(function,search-libs[,action-if-found][,action-if-not-found][,other-libraries])
# same as AC_SEARCH_LIBS, but does not try any library unless $libsearch==true
# Should be used instead of AC_SEARCH_LIBS in order to make --disable-libsearch work correctly
AC_DEFUN([PLUMED_SEARCH_LIBS], [
  if test "${libsearch}" == true ; then
    AC_SEARCH_LIBS([$1],[$2],[$3],[$4],[$5])
  else
    AC_SEARCH_LIBS([$1],[],[$3],[$4],[])
  fi
])

# PLUMED_CHECK_PACKAGE(header,function,define[,library])
# first we check if the header is present. if so, we check if the given function can be found.
# if it cannot be found and the fourth argument (library) has been supplied, we look for it
# in the library. finally, we set the define flag
AC_DEFUN([PLUMED_CHECK_PACKAGE], [
    found=ko
    $3=no
    AC_CHECK_HEADER( [$1],[
      PLUMED_SEARCH_LIBS([$2],[$4],[found=ok])
    ])
    if test $found == ok ; then
       AC_DEFINE([$3])
       $3=yes
    else
       AC_MSG_WARN([cannot enable $3]) 
    fi
])

# PLUMED_CHECK_CXX_PACKAGE(name,program,define[,library])
# Similar to PLUMED_CHECK_PACKAGE but suitable for checking C++ libraries.
# Name is just a string used in the configure log
# program is a string containing a short C++ program that should compiled.
# If present and necessary, library is also linked.
AC_DEFUN([PLUMED_CHECK_CXX_PACKAGE], [
    found=ko
    $3=no
    if test "${libsearch}" == true ; then
      testlibs="$4"
    else
      testlibs=""
    fi
    for testlib in "" $testlibs
    do
      save_LIBS="$LIBS"
      if test -n "$testlib" ; then
        AC_MSG_CHECKING([$1 with -l$testlib])
        LIBS="-l$testlib $LIBS"
      else
        AC_MSG_CHECKING([$1 without extra libs])
      fi
      AC_LINK_IFELSE([AC_LANG_SOURCE([$2])],
          [found=ok]
          [AC_MSG_RESULT([yes])],
          [AC_MSG_RESULT([no])]
      )
      if test $found == ok ; then
        break
      fi
      LIBS="$save_LIBS"
    done
    if test $found == ok ; then
      AC_DEFINE([$3])
      $3=yes
    else
      AC_MSG_WARN([cannot enable $3]) 
      LIBS="$save_LIBS"
    fi
])

##################################################################

AC_MSG_NOTICE([Optional modules are disabled by default])
rm -f src/*.on src/*.off

AC_ARG_ENABLE([modules],
  AS_HELP_STRING([--enable-modules], [all/none/reset or : separated list such as +crystallization:-bias default: reset]),
  [
   rm -f src/*.on src/*.off
   if test "${enableval}" == yes ; then
     enableval=reset
   fi
   if test "${enableval}" == no ; then
     enableval=none
   fi
   for modules_mod in `echo "${enableval}" | sed 's/:/ /g' | sed 's/+/ +/g' | sed 's/-/ -/g'`
   do
     case "$modules_mod" in
     (none)
       AC_MSG_NOTICE([Disabling all optional modules])
       rm -f src/*.off src/*.on
       touch $(grep default-on src/*/module.type | sed "s|/module.type:.*|.off|") ;;
     (all)
       AC_MSG_NOTICE([Enabling all optional modules])
       rm -f src/*.off src/*.off
       touch $(grep default-off src/*/module.type | sed "s|/module.type:.*|.on|") ;;
     (reset)
        AC_MSG_NOTICE([Resetting modules to default])
        rm -f src/*.on src/*.off ;;
     (-*)
       modules_mod=`echo "${modules_mod}" | sed "s|-||"`
       if test ! -f src/$modules_mod/module.type ; then
         AC_MSG_ERROR([trying to remove module $modules_mod which does not exist])
       fi
       AC_MSG_NOTICE([Disabling module ${modules_mod}])
       rm -f src/$modules_mod.on
       touch src/$modules_mod.off ;;
     (+*)
       modules_mod=`echo "${modules_mod}" | sed "s|+||"`
       if test ! -f src/$modules_mod/module.type ; then
         AC_MSG_ERROR([trying to add module $modules_mod which does not exist])
       fi
       AC_MSG_NOTICE([Enabling module ${modules_mod}])
       rm -f src/$modules_mod.off
       touch src/$modules_mod.on ;;
     (*)
       if test ! -f src/$modules_mod/module.type ; then
         AC_MSG_ERROR([trying to add module $modules_mod which does not exist])
       fi
       AC_MSG_NOTICE([Enabling module ${modules_mod}])
       rm -f src/$modules_mod.off
       touch src/$modules_mod.on ;;
     esac
   done
  ]
  ,
  []
)


# set enable flags for ./configure

# This can be disabled when configuring within a package manager
# such as macports to make sure that only libraries explicitly
# listed in LDFLAGS are linked
PLUMED_CONFIG_ENABLE([libsearch],[search for libraries],[yes])

# This can be disabled to avoid the extra tests for static patching.
# In the future, when this is disabled, we could disable the "--static" flag
# of "plumed patch"
PLUMED_CONFIG_ENABLE([static-patch],[allow statically linking plumed],[yes])

PLUMED_CONFIG_ENABLE([doc],[documentation],[yes])
PLUMED_CONFIG_ENABLE([pdfdoc],[pdf version of the manual],[no])
PLUMED_CONFIG_ENABLE([debug],[debugging],[no])
PLUMED_CONFIG_ENABLE([gcov],[gcov to estimate code coverage],[no])
PLUMED_CONFIG_ENABLE([basic-warnings],[basic warnings],[yes])
PLUMED_CONFIG_ENABLE([fussy],[fussy warnings],[no])
PLUMED_CONFIG_ENABLE([debug-glibcxx],[enable boundary check],[no])
PLUMED_CONFIG_ENABLE([shared],[shared libs],[yes])
PLUMED_CONFIG_ENABLE([dependency-tracking],[dependency tracking],[yes])
PLUMED_CONFIG_ENABLE([rpath],[store rpath],[no])
PLUMED_CONFIG_ENABLE([absolute-soname],[store absolute soname (Linux only - this is the default behavior on OSX). Only enable for testing!],[no])
PLUMED_CONFIG_ENABLE([absolute-install-name],[store absolute relative (OSX only - disable to have a behavior similar to Linux). Only disable for testing!],[yes])
PLUMED_CONFIG_ENABLE([loader-path],[use @loader_path to find libplumedKernel.dylib (OSX only)],[yes])
PLUMED_CONFIG_ENABLE([bsymbolic],[use -Bsymbolic flag in making shared libraries (Linux only)],[yes])
PLUMED_CONFIG_ENABLE([ld-r],[group object files],[yes])
PLUMED_CONFIG_ENABLE([ar-cr],[use ar to build libplumedWrapper.a],[yes])
PLUMED_CONFIG_ENABLE([static-archive],[try to build libplumed.a for static linking],[yes])
PLUMED_CONFIG_ENABLE([asmjit],[enable embedded asmjit],[no])
PLUMED_CONFIG_ENABLE([mpi],[search for mpi],[yes])
PLUMED_CONFIG_ENABLE([external-lapack],[search for external lapack],[yes])
PLUMED_CONFIG_ENABLE([external-blas],[search for external blas],[yes])
PLUMED_CONFIG_ENABLE([molfile-plugins],[use molfile_plugins],[yes])
PLUMED_CONFIG_ENABLE([external-molfile-plugins],[search for external molfile_plugins],[yes])
PLUMED_CONFIG_ENABLE([zlib],[search for zlib],[yes])
PLUMED_CONFIG_ENABLE([readdir-r],[search for readdir_r (threadsafe)],[no])
PLUMED_CONFIG_ENABLE([cregex],[search for C regular expression],[yes])
PLUMED_CONFIG_ENABLE([dlopen],[search for dlopen],[yes])
PLUMED_CONFIG_ENABLE([rtld_default],[search for RTLD_DEFAULT macro],[yes])
PLUMED_CONFIG_ENABLE([chdir],[search for chdir function],[yes])
PLUMED_CONFIG_ENABLE([subprocess],[search for functions needed to manage a subprocess],[yes])
PLUMED_CONFIG_ENABLE([getcwd],[search for getcwd function],[yes])
PLUMED_CONFIG_ENABLE([execinfo],[search for execinfo],[yes])
PLUMED_CONFIG_ENABLE([gsl],[search for gsl],[yes])
PLUMED_CONFIG_ENABLE([xdrfile],[search for xdrfile],[yes])
PLUMED_CONFIG_ENABLE([boost_graph],[search for boost graph],[no])
PLUMED_CONFIG_ENABLE([boost_serialization],[search for boost serialization],[no])
PLUMED_CONFIG_ENABLE([fftw],[search for fftw],[yes])
PLUMED_CONFIG_ENABLE([python],[search for python],[yes])
PLUMED_CONFIG_ENABLE([af_ocl],[search for arrayfire_ocl],[no])
PLUMED_CONFIG_ENABLE([af_cuda],[search for arrayfire_cuda],[no])
PLUMED_CONFIG_ENABLE([af_cpu],[search for arrayfire_cpu],[no])

AC_ARG_VAR(SOEXT,[extension of dynamic libraries (so/dylib)])
AC_ARG_VAR(STATIC_LIBS,[variables that should be linked statically directly to MD code - configure will add here -ldl if necessary ])
AC_ARG_VAR(LDSHARED,[command for linking shared library - configure will use CXX plus the proper flags ])
AC_ARG_VAR(PYTHON_BIN,[python executable (e.g. python2.7 or /opt/local/bin/python2.7) - default: search for a python executable])
AC_ARG_VAR(BASH_COMPLETION_DIR,[path where bash completion will be installed - default: search with pkg-config])
AC_ARG_VAR(MPIEXEC,[command to run mpi programs in tests - default not specified, which means use PLUMED_MPIRUN env var at runtime for backward compatibility])


# by default use -O flag
# we override the autoconf default (-g) because in release build we do not want to
# include symbol information (obj files are huge)
if test -z "$CXXFLAGS"
then
  CXXFLAGS=-O3
fi

# this is a list of common compilers
compilers="g++ c++ cxx icpc"
c_compilers="gcc cc icc"

# on OSX, prefer clang++
case `(uname)` in
  (Darwin)
    compilers="clang++ $compilers"
    c_compilers="clang $c_compilers"
  ;;
esac

# if searching for MPI, try first mpi-like compilers
if test $mpi == true ; then
compilers="mpic++ mpicxx mpiicpc openmpic++ openmpicxx $compilers"
fi

# do the actual search
AC_PROG_CXX([$compilers])
AC_PROG_CC([$c_compilers])

# also setup Fortran compiler
# this is optional, and can be used in the late part of this 
# script to verify that fortran can indeed link properly the
# a c++ library
AC_PROG_FC

# we use C++ for all the autoconf tests
AC_LANG(C++)

if test -z "$LDSHARED" ; then
  LDSHARED="$CXX"
fi

AC_MSG_NOTICE([Initial CXX:         $CXX])
AC_MSG_NOTICE([Initial CXXFLAGS:    $CXXFLAGS])
AC_MSG_NOTICE([Initial CPPFLAGS:    $CPPFLAGS])
AC_MSG_NOTICE([Initial CFLAGS:      $CFLAGS])
AC_MSG_NOTICE([Initial LDFLAGS:     $LDFLAGS])
AC_MSG_NOTICE([Initial LIBS:        $LIBS])
AC_MSG_NOTICE([Initial STATIC_LIBS: $STATIC_LIBS])
AC_MSG_NOTICE([Initial LD:          $LD])
AC_MSG_NOTICE([Initial LDSHARED:    $LDSHARED])
AC_MSG_NOTICE([Initial SOEXT:       $SOEXT])

# check C++ flags
if test $shared == true
then
  PLUMED_CHECK_CXXFLAG([-fPIC])
  PLUMED_CHECK_CFLAG([-fPIC])
fi

if test $basic_warnings == true
then
  PLUMED_CHECK_CXXFLAG([-Wall])
  PLUMED_CHECK_CXXFLAG([-pedantic])
  PLUMED_CHECK_CXXFLAG([-std=c++11])
fi

if test $debug == true
then
  PLUMED_CHECK_CXXFLAG([-g])
fi

if test $gcov == true
then
  CXX="$CXX --coverage"
  LDSHARED="$LDSHARED --coverage"
  CFLAGS="$CFLAGS --coverage"
  PLUMED_CHECK_LDFLAGS([--coverage])
  STATIC_LIBS="$STATIC_LIBS --coverage"
fi

if test $fussy == true
then
  PLUMED_CHECK_CXXFLAG([-Wextra])
  PLUMED_CHECK_CXXFLAG([-Wfloat-equal])
  PLUMED_CHECK_CXXFLAG([-Wwrite-strings])
  PLUMED_CHECK_CXXFLAG([-Wpointer-arith])
  PLUMED_CHECK_CXXFLAG([-Wcast-qual])
  PLUMED_CHECK_CXXFLAG([-Wcast-align])
  PLUMED_CHECK_CXXFLAG([-Wconversion])
  PLUMED_CHECK_CXXFLAG([-Wredundant-delcs])
  PLUMED_CHECK_CXXFLAG([-Wvariadic-macros])
  PLUMED_CHECK_CXXFLAG([-Wold-style-cast])
fi

AC_MSG_CHECKING([whether $CXX declares c++11 support])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
#if __cplusplus <= 199711L
this_compiler_does_not_support_cxx11
#endif
])],
  [AC_MSG_RESULT([yes])];
  support_cxx11=true,
  [AC_MSG_RESULT([no]) ;
  support_cxx11=false])

if test "$support_cxx11" = false
then
   AC_MSG_WARN([C++11 support is required as of PLUMED 2.4])
   AC_MSG_WARN([Your compiler appears not to support C++11])
   AC_MSG_WARN([Please change compiler or make sure that everything works correctly])
fi

AC_MSG_CHECKING([whether C++ library supports C++11 exceptions])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
#include <new>
#include <functional>
#include <memory>
#include <system_error>
#include <ios>
#include <future>
// capable to catch:
void func(void) {
  try{
  } catch(std::bad_array_new_length &) {
  } catch(std::bad_function_call &) {
  } catch(std::bad_weak_ptr &) {
  } catch(std::ios_base::failure &e) {
    if(e.code().category()==std::generic_category()) {}
    else if(e.code().category()==std::system_category()) {}
    else if(e.code().category()==std::iostream_category()) {}
    else if(e.code().category()==std::future_category()) {}
  } catch(std::system_error &e) {
    if(e.code().category()==std::generic_category()) {}
    else if(e.code().category()==std::system_category()) {}
    else if(e.code().category()==std::iostream_category()) {}
    else if(e.code().category()==std::future_category()) {}
  }
// capable to throw:
  auto a=std::bad_array_new_length();
  auto b=std::bad_function_call();
  auto c=std::bad_weak_ptr();
  auto d=std::system_error(10,std::generic_category(),"msg");
  auto e=std::system_error(10,std::system_category(),"msg");
  auto f=std::system_error(10,std::iostream_category(),"msg");
  auto g=std::system_error(10,std::future_category(),"msg");
  auto h=std::ios_base::failure("msg",std::error_code(10,std::generic_category()));
  auto i=std::ios_base::failure("msg",std::error_code(10,std::system_category()));
  auto j=std::ios_base::failure("msg",std::error_code(10,std::iostream_category()));
  auto k=std::ios_base::failure("msg",std::error_code(10,std::future_category()));
}
])],
[
  AC_MSG_RESULT([yes]);
  AC_DEFINE([__PLUMED_LIBCXX11])
],
[AC_MSG_RESULT([no])
])

AC_SUBST(disable_dependency_tracking)

if test "$dependency_tracking" = true
then
  AC_MSG_CHECKING([whether $CXX can generate dependency file with -MM -MF])
  dependency=ko
  echo "#include \"conftest1.h\"" > conftest.cpp
  echo "#include \"conftest2.h\"" > conftest1.h
  echo "/* */" > conftest2.h
  $CXX $CXXFLAGS -c -MM -MFconftest.d conftest.cpp 1> /dev/null 2> /dev/null
  grep conftest2 conftest.d 1> /dev/null 2>/dev/null && dependency=ok
  if test "$dependency" = ok ; then
    AC_MSG_RESULT([yes])
    disable_dependency_tracking=no
  else
    AC_MSG_RESULT([no])
    disable_dependency_tracking=yes
  fi
else
  disable_dependency_tracking=yes
fi

if test "$disable_dependency_tracking" = yes ; then
  AC_MSG_WARN([dependencies tracking disabled - always make clean before make])
else
  AC_MSG_NOTICE([dependency tracking enabled])
fi




#### Compulsory libraries ####
# some of them might be made optional if we find that are not available in some system
AC_MSG_NOTICE([Now we will check compulsory headers and libraries])
AC_CHECK_HEADER([dirent.h],     [ ], [AC_MSG_ERROR([compulsory header not found])] )
AC_CHECK_FUNC(  [readdir],      [ ], [AC_MSG_ERROR([compulsory function not found])] )

# Then check for blas. This is a bit complicated because we want to allow
# either the version with underscore or the one without
blas_found=
lapack_found=

# external lapack can only work with external blas
# thus, if external blas are disabled also external lapack should be disabled
if test "$external_blas" == false && test "$external_lapack" == true ; then
  AC_MSG_NOTICE([Internal blas can only be used with internal lapack])
  AC_MSG_NOTICE([Will not search for external lapack])
  external_lapack=false
fi

# first look for blas
if test "$external_blas" == true ; then
  AC_MSG_CHECKING([whether dgemv can be linked with no library])
  AC_LINK_IFELSE([AC_LANG_CALL([], [dgemv])],[
    AC_MSG_RESULT([yes])
    blas_found=nounderscore
    ],[
    AC_MSG_RESULT([no])
    AC_MSG_CHECKING([whether dgemv_ can be linked with no library])
    AC_LINK_IFELSE([AC_LANG_CALL([], [dgemv_])],[
      AC_MSG_RESULT([yes])
      blas_found=underscore
    ],[
      AC_MSG_RESULT([no])
      PLUMED_SEARCH_LIBS([dgemv],[blas],[blas_found=nounderscore],[
        PLUMED_SEARCH_LIBS([dgemv_],[blas],[blas_found=underscore])
      ])
    ])
  ])
fi

# if not found, then use internal lapack and blas
if test -z "$blas_found" ; then
AC_MSG_WARN([using internal lapack and blas, could be inefficient])
fi

# if found, also look for external lapack
if test -n "$blas_found" ; then

AC_DEFINE([__PLUMED_HAS_EXTERNAL_BLAS])

if test "$external_lapack" == true ; then
# Then we look for lapack using same underscoring
case "$blas_found" in
(underscore) search_for=dsyevr_ ;;
(nounderscore) search_for=dsyevr ;;
esac
PLUMED_SEARCH_LIBS( [$search_for],[lapack],[lapack_found=yes])
fi

# if not found, then use internal lapack with external blas
if test -z "$lapack_found" ; then
AC_MSG_WARN([using internal lapack, could be inefficient])
else
AC_DEFINE([__PLUMED_HAS_EXTERNAL_LAPACK])
fi

fi 

# in case external blas have been found, take note of their underscoring
# notice that this applies also when external blas are used with internal lapack
# in the latter case, also (internal) lapack names will be underscored consistently
if test "$blas_found" == nounderscore
then
  AC_DEFINE([F77_NO_UNDERSCORE])
fi

if test -n "$blas_found" ; then

AC_DEFUN([PLUMED_TEST_SDOT],[
#if ! defined(F77_NO_UNDERSCORE)
#define sdot sdot_
#endif
extern "C"{
$1 sdot(int *n, float *dx, int *incx, float *dy, int *incy);
}

int main(){
  int size=5;
// notice that [[]] is required to cheat autoconf
  float af[[5]],bf[[5]];
  for(unsigned i=0;i<size;i++){
    af[[i]]=size;
    bf[[i]]=size-i;
  }
  int inca=1;
  int incb=1;
  $1 f=sdot(&size,af,&inca,bf,&incb)-75;
  if(f*f<1e-10) return 0;
  else return 1;
}
])


AC_MSG_CHECKING([whether sdot returns float])
AC_RUN_IFELSE([AC_LANG_SOURCE([PLUMED_TEST_SDOT([float])])],
  [ sdot_returns_float=yes ; AC_MSG_RESULT([yes]) ],
  [ sdot_returns_float=no ;  AC_MSG_RESULT([no]) ],
  [ AC_MSG_RESULT([not checked (cross compiling)]) ]
)

if test "$sdot_returns_float" == no ;
then

AC_MSG_CHECKING([whether sdot returns double])
AC_RUN_IFELSE([AC_LANG_SOURCE([PLUMED_TEST_SDOT([double])])],
  [ sdot_returns_double=yes ; AC_MSG_RESULT([yes]) ],
  [ sdot_returns_double=no ;  AC_MSG_RESULT([no]) ],
  [ AC_MSG_RESULT([not checked (cross compiling)]) ]
)

fi

if test "$sdot_returns_double" = yes ; then
  AC_MSG_NOTICE([Setting workaround for blas float functions returning double])
  AC_DEFINE(__PLUMED_BLAS_RETURNS_FLOAT,double)
  if test -n "$lapack_found" ; then
    AC_MSG_NOTICE([Setting workaround for lapack float functions returning double])
    AC_DEFINE(__PLUMED_LAPACK_RETURNS_FLOAT,double)
  fi
else if test "$sdot_returns_float" = no && test "$sdot_returns_double" = no ;
then
  AC_MSG_WARN([There is a problem with your blas implementation])
fi

fi

fi

#### End of compulsory libraries ####

#### Optional libraries ####
AC_MSG_NOTICE([Now we will check for optional headers and libraries])

#############################################################
# I add the possibility to completely remove molfile_plugins
# I would like to be 100% that the molfile module compiles
# correctly on all machines
# In case of problem, it is sufficient to configure with
# ./configure --disable-molfile-plugins
# GB
#############################################################

if test $molfile_plugins == true ; then

# Check for molfile_plugins and use internal fallback if not found. TG

# We always have molfile, now
AC_DEFINE([__PLUMED_HAS_MOLFILE_PLUGINS])

  if test "$external_molfile_plugins" == true ; then
    PLUMED_CHECK_PACKAGE([libmolfile_plugin.h],[molfile_dcdplugin_init],[__PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS],molfile_plugin)
    if test "$__PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS" != yes ; then
      AC_MSG_WARN([using internal molfile_plugins, which only support dcd/xtc/trr/trj/crd files])
    else
      AC_DEFINE([__PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS])
    fi
  fi
fi

# this is special and is also attached to STATIC_LIBS
# this flag should be used also when linking MD engines to allow plumed
# to be loaded later
AC_CHECK_LIB([dl],dlopen, [STATIC_LIBS="-ldl $STATIC_LIBS"] [LIBS="-ldl $LIBS"])

mpi_found=ko
# optional libraries follow
if test $mpi == true ; then
  PLUMED_CHECK_PACKAGE([mpi.h],[MPI_Init],[__PLUMED_HAS_MPI])
  if test "$__PLUMED_HAS_MPI" = yes; then
    mpi_found=ok
  fi
else 
  mpi_found=ko
fi

# search for openmp is automatically disabled by autoconf
# when configuring with --disable-openmp
AC_OPENMP

if test $asmjit == true ; then
  found=ko
# asmjit calls clock_gettime and thus should be linked to rt on Linux
  if test `(uname)` = Linux ; then
    PLUMED_SEARCH_LIBS([clock_gettime],[rt],[found=ok])
  else
    found=ok
  fi
  if test "$found" = ok ; then
    AC_MSG_NOTICE([Enabling embedded asmjit])
    AC_DEFINE([__PLUMED_HAS_ASMJIT])
  else
    AC_MSG_WARN([cannot link clock_gettime on this Linux, asmjit will not be enabled])
  fi
fi


if test $readdir_r == true ; then
  PLUMED_CHECK_PACKAGE([dirent.h],[readdir_r],[__PLUMED_HAS_READDIR_R])
fi
if test $cregex == true ; then
  PLUMED_CHECK_PACKAGE([regex.h],[regcomp],[__PLUMED_HAS_CREGEX])
fi
if test $dlopen == true ; then
  PLUMED_CHECK_PACKAGE([dlfcn.h],[dlopen],[__PLUMED_HAS_DLOPEN])
fi

if test $rtld_default == true ; then
  PLUMED_CHECK_CXX_PACKAGE([RTLD_DEFAULT],[
#include <dlfcn.h>
int
main ()
{
  void* f=dlsym(RTLD_DEFAULT,"path");
  return 0;
}
  ], [__PLUMED_HAS_RTLD_DEFAULT])
fi

if test $chdir == true ; then
  PLUMED_CHECK_PACKAGE([unistd.h],[chdir],[__PLUMED_HAS_CHDIR])
fi

if test $subprocess == true ; then
  PLUMED_CHECK_CXX_PACKAGE([subprocess],[
/* random program calling all the functions needed to manage a subprocess */
#include <cstdio>
#include <unistd.h>
#include <csignal>

int
main ()
{
// notice that [[]] is required to cheat autoconf
  int cp[[2]];
  FILE* fp;
  char* arr[[3]];
  arr[[0]]=NULL;
  arr[[1]]=NULL;
  arr[[2]]=NULL;
  arr[[3]]=NULL;
  if(pipe(cp)>=0)
  if(fork()>=0)
  if(close(1)>=0)
  if(dup(cp[[1]])>=0) {
    fp=fdopen(cp[[0]],"w");
    execv(arr[[0]],arr);
  }
  auto p=fork();
  kill(p,SIGINT);
  kill(p,SIGCONT);
  kill(p,SIGSTOP);
  return 0;
}

  ], [__PLUMED_HAS_SUBPROCESS])
fi

if test $getcwd == true ; then
  PLUMED_CHECK_PACKAGE([unistd.h],[getcwd],[__PLUMED_HAS_GETCWD])
fi

if test $execinfo == true ; then
  PLUMED_CHECK_PACKAGE([execinfo.h],[backtrace],[__PLUMED_HAS_EXECINFO])
fi
if test $zlib == true ; then
  PLUMED_CHECK_PACKAGE([zlib.h],[gzopen],[__PLUMED_HAS_ZLIB],[z])
fi

if test $gsl == true ; then
  found=ko
  PLUMED_SEARCH_LIBS([cblas_dgemv],[gslcblas], [
    AC_CHECK_HEADER(  [gsl/gsl_vector.h], [
      PLUMED_SEARCH_LIBS([gsl_vector_alloc],[gsl],[found=ok])
    ])
  ])
  if test $found == ok ; then
    AC_DEFINE([__PLUMED_HAS_GSL])
  else
    AC_MSG_WARN([cannot enable __PLUMED_HAS_GSL])
  fi
fi

if test $xdrfile == true ; then
  PLUMED_CHECK_PACKAGE([xdrfile/xdrfile_xtc.h],[write_xtc],[__PLUMED_HAS_XDRFILE],[xdrfile])
fi

if test $boost_graph == true ; then
  PLUMED_CHECK_CXX_PACKAGE([boost graph],[
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_matrix.hpp>
int
main ()
{
  boost::adjacency_matrix<boost::directedS> a(1);
  ;
  return 0;
}
  ], [__PLUMED_HAS_BOOST_GRAPH])
fi

if test $boost_serialization == true ; then
  PLUMED_CHECK_CXX_PACKAGE([boost serialization],[
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
int main() {
    std::ofstream ofs("filename");
    boost::archive::text_oarchive oa(ofs);
    return 0;
}
  ], [__PLUMED_HAS_BOOST_SERIALIZATION],[boost_serialization boost_serialization-mt])
# notice: macports install libraries with -mt suffix
fi

if test $fftw == true ; then
  PLUMED_CHECK_PACKAGE([fftw3.h],[fftw_execute],[__PLUMED_HAS_FFTW],[fftw3])
fi

if test $python == true  ; then
# if PYTHON_BIN is defined, it is expected to be the full path to python
# Otherwise, search from a list of names:
  if test -z "$PYTHON_BIN" ; then
    AC_CHECK_PROGS([PYTHON_BIN],[python])
  fi
  if test -n "$PYTHON_BIN"
  then
    AC_MSG_NOTICE([Python executable is $PYTHON_BIN])
    AC_MSG_CHECKING([support for required python modules (python3, setuptools, cython, subprocess, os, shutil)])
testimport="
from setuptools import setup
from setuptools import Extension
import subprocess
import os
import os.path
import sys
from shutil import copyfile
from Cython.Build import cythonize
if sys.version_info < (3,):
    raise ImportError('PLUMED 2.6 only supports Python 3')
"
    if echo "$testimport" | "$PYTHON_BIN" 1>/dev/null 2>/dev/null;  then
       AC_MSG_RESULT([yes])
       AC_DEFINE([__PLUMED_HAS_PYTHON])
    else
       AC_MSG_RESULT([no])
       AC_MSG_WARN([cannot enable python interface])
       PYTHON_BIN=
    fi
  else
    AC_MSG_WARN([cannot enable python interface])
  fi
fi

if test "$af_ocl" == true ; then
  PLUMED_CHECK_PACKAGE([arrayfire.h],[af_is_double],[__PLUMED_HAS_ARRAYFIRE],[afopencl])
fi
if test "$af_cuda" == true ; then
  PLUMED_CHECK_PACKAGE([arrayfire.h],[af_is_double],[__PLUMED_HAS_ARRAYFIRE],[afcuda])
fi
if test "$af_cpu" == true ; then
  PLUMED_CHECK_PACKAGE([arrayfire.h],[af_is_double],[__PLUMED_HAS_ARRAYFIRE],[afcpu])
fi

# in non-debug mode, add -DNDEBUG
if test "$debug" == false ; then
  AC_MSG_NOTICE([Release mode, adding -DNDEBUG])
  AC_DEFINE([NDEBUG])
fi

# in debug-glibcxx mode, add -D_GLIBCXX_DEBUG
if test "$debug_glibcxx" == true ; then
  AC_MSG_NOTICE([Check boundaries, adding -D_GLIBCXX_DEBUG])
  AC_DEFINE([_GLIBCXX_DEBUG])
fi

# this is necessary in many MPI implementations
# I leave it by default, since it seems harmless
AC_DEFINE([_REENTRANT])

#### Options for dynamic library to work properly ####
AC_SUBST(SOEXT)
AC_SUBST(LDSHARED)
# these are libraries that should be linked also to MD engines
AC_SUBST(STATIC_LIBS)
# python executable
AC_SUBST(PYTHON_BIN)

AC_SUBST(MPIEXEC)

if test "$shared" == true ; then
  case `(uname)` in
  (Darwin)
    AC_MSG_NOTICE([*** Special settings for dynamic libraries on OSX ***])
    AC_MSG_NOTICE([Dynamic library extension is 'dylib'])
    AC_MSG_NOTICE([LDSHARED needs special flags])
    SOEXT=dylib
    LDSHARED="$LDSHARED -dynamiclib -Wl,-headerpad_max_install_names"
    if test "$rpath" = true ; then
      AC_MSG_NOTICE([Switching off rpath on OSX])
      rpath=false
    fi
  ;;
  (Linux)
    AC_MSG_NOTICE([*** Special settings for dynamic libraries on Linux ***])
    AC_MSG_NOTICE([Dynamic library extension is 'so'])
    AC_MSG_NOTICE([LDSHARED and LDFLAGS need special flags])
    SOEXT=so
    LDSHARED="$LDSHARED -shared"
    PLUMED_CHECK_LDFLAGS([-rdynamic])
    if test "$bsymbolic" == true ; then
      PLUMED_CHECK_LDFLAGS([-Wl,-Bsymbolic])
    fi
  ;;
  (*)
    AC_MSG_NOTICE([*** Dynamic library only enabled on OSX and Linux ***])
  esac
fi

# check linking of runtime library
if test -n "$SOEXT"
then
  AC_MSG_NOTICE([Using LDSHARED='$LDSHARED'])
  AC_MSG_NOTICE([Using LDFLAGS='$LDFLAGS'])
  AC_MSG_CHECKING([whether LDSHARED can create dynamic libraries])
  rm -f conftest.*
  echo "void f(void){ return;}" > conftest.cpp
  $CXX $CXXFLAGS $CPPFLAGS -c conftest.cpp 1>/dev/null 2>/dev/null
  $LDSHARED $LDFLAGS conftest.o -o conftest.$SOEXT 1>/dev/null 2>/dev/null
  if test -f conftest.$SOEXT
  then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    AC_MSG_WARN([dynamic library will be disabled])
    SOEXT=
  fi
  if test -n "$SOEXT" && test "$rpath" = true ; then
    rm -f conftest.$SOEXT
    readelf=""
    AC_CHECK_PROG([readelf],[readelf],[found])
    if test "$readelf" == found ; then
      test_LDSHARED="$LDSHARED -Wl,-R -Wl,/some/random/dir/"
      $test_LDSHARED $LDFLAGS conftest.o -o conftest.$SOEXT 1>/dev/null 2>/dev/null
      if readelf -d conftest.$SOEXT | grep RPATH | grep -q /some/random/dir ; then
        LDSHARED="$LDSHARED -Wl,-R -Wl,\"$LIBRARY_PATH\""
      else
        AC_MSG_WARN([-R option seems not working, disabling rpath])
      fi
    else
      AC_MSG_WARN([readelf not available, no way to set rpath])
    fi 
  fi
  rm -f conftest.*
fi
#### Options for dynamic library to work properly ####

AC_SUBST(make_doc)
make_doc=no
if test "$doc" == true
then

make_doc=yes

### Look for doxygen
AC_CHECK_PROG([doxygen],[doxygen],[found])
if test "$doxygen" == found
then
  doxygen_version=`doxygen --version | awk 'BEGIN{FS="."}{if($1>1 || ($1==1 && $2>=8)) print "ok"}'`
  if test "$doxygen_version" == ok
  then
    AC_MSG_NOTICE([Doxygen version is fine])
  else
    AC_MSG_WARN([Doxygen version is <1.8])
    make_doc=no
  fi
  AC_CHECK_PROG([dot],[dot],[found])
  if test "$dot" != found
  then
    AC_MSG_WARN([You will not be able to see diagrams in the manual])
  fi
else
  make_doc=no
fi
fi

if test "$make_doc" = yes
then
  AC_MSG_NOTICE([Manuals will be generated])
else
  AC_MSG_NOTICE([Manuals will not be generated])
fi

AC_SUBST(make_pdfdoc)
make_pdfdoc=""
if test "$pdfdoc" == true && test "$make_doc" == yes
then
  AC_MSG_NOTICE([A PDF version of the manual will be generated])
  make_pdfdoc=yes
else
  AC_MSG_NOTICE([A PDF version of the manual will not be generated])
  make_pdfdoc=no
fi

### Look for xxd
AC_CHECK_PROG([xxd],[xxd],[found])
if test "$xxd" != found
then
  AC_MSG_ERROR([xxd should be installed for PLUMED to compile properly])
fi

AC_SUBST(program_can_run)
program_can_run=""
AC_MSG_CHECKING([whether a program can be run on this machine])
AC_RUN_IFELSE([AC_LANG_SOURCE([
#ifdef __PLUMED_HAS_MPI
#include <mpi.h>
#endif
// notice that [[]] is required to cheat autoconf
int main(int argc,char*argv[[]]){
#ifdef __PLUMED_HAS_MPI
// this emulates what happens when plumed
// is compiled with mpi and invoked with --no-mpi
  if(argc==10){
    MPI_Init(&argc,&argv);
    return MPI_Finalize();
  }
#endif
  return 0;
}
])],
  [ program_can_run=yes ; AC_MSG_RESULT([yes]) ],
  [ program_can_run=no ; AC_MSG_RESULT([no]) ],
  [ program_can_run=no ; AC_MSG_RESULT([no (cross compiling)]) ]
)

if test $mpi_found == ok ; then
AC_SUBST(program_can_run_mpi)
program_can_run_mpi=""
AC_MSG_CHECKING([whether a program compiled with mpi can be run on this machine])
AC_RUN_IFELSE([AC_LANG_SOURCE([
#ifdef __PLUMED_HAS_MPI
#include <mpi.h>
#endif
// notice that [[]] is required to cheat autoconf
int main(int argc,char*argv[[]]){
#ifdef __PLUMED_HAS_MPI
  MPI_Init(&argc,&argv);
  return MPI_Finalize();
#endif
  return 0;
}
])],
  [ program_can_run_mpi=yes ; AC_MSG_RESULT([yes]) ],
  [ program_can_run_mpi=no ; AC_MSG_RESULT([no]) ],
  [ program_can_run_mpi=no ; AC_MSG_RESULT([no (cross compiling)]) ]
)
fi

if test $mpi_found == ok ; then
  if test -n "$MPIEXEC" ; then
    AC_MSG_NOTICE([Regtest suite will use $MPIEXEC command to run MPI tests])
  else
    AC_MSG_NOTICE([Regtest suite will use env var PLUMED_MPIRUN to run MPI tests (default: mpirun)])
  fi
fi

if test "$SOEXT" == "dylib" ; then
  use_absolute_soname=yes
else
  use_absolute_soname=no
fi

AC_SUBST(use_absolute_soname)
if test "$absolute_soname" == true ; then
  if test "$SOEXT" == "dylib" ; then
    AC_MSG_WARN([--enable-absolute-soname has no effect on OSX])
  else
    AC_MSG_NOTICE([enabling absolute soname. Full path will be hardcoded in plumed library soname])
    use_absolute_soname=yes
  fi
fi

if test "$absolute_install_name" == false ; then
  if test "$SOEXT" == "so" ; then
    AC_MSG_WARN([--disable-absolute-install-name has no effect on Linux])
  else
    AC_MSG_NOTICE([enabling relative install_name. You will have to set DYLD_LIBRARY_PATH yor plumed libraries to be found at runtime])
    use_absolute_soname=no
  fi
fi

AC_SUBST(use_loader_path)
if test "$SOEXT" == "dylib" ; then
  if test "$loader_path" == true ; then
    use_loader_path=yes
  else
    use_loader_path=no
  fi
else
  use_loader_path=no
fi

#### This further tests are required to allow linking with non c++ compiler
AC_MSG_NOTICE([PLUMED seems to be configured properly!])
AC_MSG_NOTICE([**************************])
AC_SUBST(LD_RO)
LD_RO=
if test "$ld_r" == true ; then
for test_LD_RO in "$($CXX --print-prog-name=ld) -r -o" "$CXX -Wl,-r -o" "ld -r -o"
do
  AC_MSG_CHECKING([whether C++ objects can be grouped with $test_LD_RO])
  
  rm -f conftest-*
  
  cat << EOF > conftest-main.cpp
  void f(void);
  int main(int argc,char**argv){ f(); return 0; }
EOF
  cat << EOF > conftest-f.cpp
  void g(void);
  void f(void){ g(); }
EOF
  cat << EOF > conftest-g.cpp
  void g(void){ return; }
EOF
  
  $CXX $CXXFLAGS -c conftest-main.cpp 1> /dev/null 2> /dev/null
  $CXX $CXXFLAGS -c conftest-f.cpp 1> /dev/null 2> /dev/null
  $CXX $CXXFLAGS -c conftest-g.cpp 1> /dev/null 2> /dev/null
  
  $test_LD_RO conftest-both.o conftest-f.o conftest-g.o 1> /dev/null 2> /dev/null
  
  $CXX $CXXFLAGS -o conftest.exe conftest-main.o conftest-both.o 1> /dev/null 2> /dev/null
  
  if test -f conftest.exe
  then
    AC_MSG_RESULT([yes])
    LD_RO="$test_LD_RO"
    break
  else
    AC_MSG_RESULT([no])
  fi
done

fi

AC_SUBST(AR_CR)
AR_CR=
if test "$ar_cr" == true ; then
for test_AR_CR in "$($CXX --print-prog-name=ar) cr" "ar cr"
do
  AC_MSG_CHECKING([whether static libraries can be created with $test_AR_CR])
  
  rm -f conftest-*
  
  cat << EOF > conftest-main.cpp
  void f(void);
  int main(int argc,char**argv){ f(); return 0; }
EOF
  cat << EOF > conftest-f.cpp
  void g(void);
  void f(void){ g(); }
EOF
  cat << EOF > conftest-g.cpp
  void g(void){ return; }
EOF
  
  $CXX $CXXFLAGS -c conftest-main.cpp 1> /dev/null 2> /dev/null
  $CXX $CXXFLAGS -c conftest-f.cpp 1> /dev/null 2> /dev/null
  $CXX $CXXFLAGS -c conftest-g.cpp 1> /dev/null 2> /dev/null
  
  $test_AR_CR conftest-both.a conftest-f.o conftest-g.o 1> /dev/null 2> /dev/null
  
  $CXX $CXXFLAGS -o conftest.exe conftest-main.o conftest-both.a 1> /dev/null 2> /dev/null
  
  if test -f conftest.exe
  then
    AC_MSG_RESULT([yes])
    AR_CR="$test_AR_CR"
    break
  else
    AC_MSG_RESULT([no])
  fi
done
 
fi

make_static_archive=no
AC_SUBST(make_static_archive)

if test "${static_archive}" == true ; then
  if test -z "$LD_RO" || test -z "$AR_CR" ; then
    AC_MSG_WARN([no way to create a static archive if ld -ro or ar cr do not work])
    static_archive=false
  fi
fi

if test "${static_archive}" == true ; then
  AC_MSG_CHECKING([whether static-object constructors can be linked from a static archive])

  magic_token=c1bc476d093a3a5c67b4530e6c54c633593fe9aa
  rm -f conftest-*
  
  cat << EOF > conftest-main.cpp
  void f(void);
  int main(int argc,char**argv){ f(); return 0; }
EOF
  cat << EOF > conftest-f.cpp
  void f(void){ return; }
EOF
  cat << EOF > conftest-g.cpp
#include <iostream>
  class g{
    public:
    g(){ std::cout<<"$magic_token\n"; }
  } init;
EOF

  $CXX $CXXFLAGS -c conftest-main.cpp 1> /dev/null 2> /dev/null
  $CXX $CXXFLAGS -c conftest-f.cpp 1> /dev/null 2> /dev/null
  $CXX $CXXFLAGS -c conftest-g.cpp 1> /dev/null 2> /dev/null
  
  $LD_RO conftest-both.o conftest-f.o conftest-g.o 1> /dev/null 2> /dev/null
# linking the previously merged objects should work:
  $AR_CR conftest-both.a conftest-both.o 1> /dev/null 2> /dev/null
# something like the following, instead, should not work:
#  $AR_CR conftest-both.a conftest-f.o conftest-g.o 1> /dev/null 2> /dev/null
#
  $CXX $CXXFLAGS -o conftest.exe conftest-main.o conftest-both.a 1> /dev/null 2> /dev/null
  if grep -q $magic_token ./conftest.exe ; then
    AC_MSG_RESULT([yes])
    make_static_archive=yes
  else
    AC_MSG_RESULT([no])
  fi
fi

if test "${static_patch}" == true ; then

AC_MSG_NOTICE([I will now check if C++ objects can be linked by C/Fortran compilers])
AC_MSG_NOTICE([This is relevant if you want to use plumed patch --static on a non-C++ code])

for compiler in CC FC
do
  rm -f conftest.* conftest-main.*
  eval compexe=\$$compiler
  if test -n "$compexe" ; then
    case $compiler in
    (CC)
      name=C
      cat << EOF > conftest-main.c
int main(int argc,char**argv){
  return 0;
}
EOF
      $CC -c conftest-main.c
    ;;
    (FC)
      name=FORTRAN
      cat << EOF > conftest-main.f90
       program main
       integer i
       end program main
EOF
      $FC -c conftest-main.f90
    ;;
    esac
    cat << EOF > conftest.cpp
#include <iostream>
void f(void){
  std::cout<<"ciao";return;
}
EOF
    $CXX $CXXFLAGS -c conftest.cpp 
# start search:
    found=
    if test "${libsearch}" == true ; then
      testlibs="-lstdc++ -lc++ -lmpi_cxx"
    else
      testlibs=""
    fi
    for testlib in "" $testlibs ; do
      comment=
      test -n "$testlib" && comment=" with library $testlib"
      AC_MSG_CHECKING([whether $name can link a C++ object$comment])
      $compexe conftest.o conftest-main.o $LDFLAGS $testlib $LIBS -o conftest.exe 1>/dev/null 2>/dev/null
      if test -f conftest.exe
      then
        found=yes
        AC_MSG_RESULT([yes])
        LIBS="$testlib $LIBS"
        break
      else
        AC_MSG_RESULT([no])
      fi
    done
    if test -z "$found" ; then
      AC_MSG_WARN([You might have problems linking $name programs.])
      AC_MSG_WARN([Please add c++ library to LIBS])
    fi
  else
    AC_MSG_NOTICE([$compiler compiler not configured])
  fi
  rm -f conftest.* conftest-main.*
done

else
AC_MSG_NOTICE([Static patching is disabled, thus tests required for static patching will be skipped])
fi

if test "$prefix" == NONE
then
  prefix=/usr/local
fi

pkgconfig_bin=""

AC_PATH_PROGS(pkgconfig_bin,pkg-config)

if test -z "$BASH_COMPLETION_DIR" && test -n "$pkgconfig_bin"; then
  if test "$prefix" == "$(pkg-config --variable=prefix bash-completion 2>/dev/null)"
  then
    AC_MSG_NOTICE([bash-completion is installed on the same prefix where plumed will be installed])
    if pkg-config --variable=completionsdir bash-completion 2>/dev/null >/dev/null ; then
      BASH_COMPLETION_DIR="$(pkg-config --variable=completionsdir  bash-completion 2>/dev/null)"
    fi
  else
    AC_MSG_NOTICE([bash-completion is not installed or it is installed on a different prefix])
    pkgconfig_bin=""
  fi
fi

if test -n "$BASH_COMPLETION_DIR"
then
  AC_MSG_NOTICE([bash completion for plumed will be installed in $BASH_COMPLETION_DIR])
else
  AC_MSG_NOTICE([bash completion for plumed will not be installed])
fi

AC_SUBST(BASH_COMPLETION_DIR)

AC_ARG_PROGRAM


# version modified to work in shell script instead of makefile:
program_transform_name_sh=$(echo "${program_transform_name}" | sed 's:\$\$:$:g')
program_name=$(echo plumed | sed "$program_transform_name_sh")
AC_SUBST(program_name)

if test "$(echo "$program_name" | tr '[A-Z]' '[a-z]')" != "$(echo "$program_name" | tr '[A-Z]' '[a-z]' | sed 's/wrapper$//')" ; then
  AC_MSG_ERROR([$program_name is not a valid program name (should not terminate with Wrapper)])
fi
if test "$(echo "$program_name" | tr '[A-Z]' '[a-z]')" != "$(echo "$program_name" | tr '[A-Z]' '[a-z]' | sed 's/kernel$//')" ; then
  AC_MSG_ERROR([$program_name is not a valid program name (should not terminate with Kernel)])
fi
if test "$(echo "$program_name" | tr '[A-Z]' '[a-z]')" != "$(echo "$program_name" | tr '[A-Z]' '[a-z]' | sed 's/-patch$//')" ; then
  AC_MSG_ERROR([$program_name is not a valid program name (should not terminate with -patch)])
fi
if test "$(echo "$program_name" | tr '[A-Z]' '[a-z]')" != "$(echo "$program_name" | tr '[A-Z]' '[a-z]' | sed 's/-config$//')" ; then
  AC_MSG_ERROR([$program_name is not a valid program name (should not terminate with -config)])
fi

AC_MSG_NOTICE([**** PLUMED will be installed using the following paths:])
AC_MSG_NOTICE([**** prefix: $prefix])
if test "$exec_prefix" = NONE ; then
  exec_prefix_='${prefix}'
else
  exec_prefix_="${exec_prefix}"
fi
AC_MSG_NOTICE([**** exec_prefix: $exec_prefix_])
AC_MSG_NOTICE([**** bindir: $bindir])
AC_MSG_NOTICE([**** libdir: $libdir])
AC_MSG_NOTICE([**** includedir: $includedir])
AC_MSG_NOTICE([**** datarootdir: $datarootdir])
AC_MSG_NOTICE([**** datadir: $datadir])
AC_MSG_NOTICE([**** docdir: ${datarootdir}/doc/$program_name])
AC_MSG_NOTICE([**** htmldir: $htmldir])
AC_MSG_NOTICE([**** Executable will be named $program_name])
AC_MSG_NOTICE([**** As of PLUMED 2.5, you cannot change paths anymore during "make install"])
AC_MSG_NOTICE([**** Please configure and make clean to change the prefix])

if test -z "$BASH_COMPLETION_DIR" ; then
  AC_MSG_WARN([**** Bash completion for plumed will not be installed, please add the following two lines to your bashrc])
  AC_MSG_WARN([**** _$program_name() { eval "\$($program_name --no-mpi completion 2>/dev/null)";}])
  AC_MSG_WARN([**** complete -F _$program_name -o default $program_name])
fi

if test $mpi == true; then
  if test $mpi_found == ok; then
    AC_MSG_NOTICE([**** PLUMED will be compiled using MPI])
  else 
    AC_MSG_WARN([**** PLUMED will NOT be compiled using MPI because MPI have not been found!])
  fi
else
    AC_MSG_NOTICE([**** PLUMED will be compiled without MPI])
fi

if test $program_can_run == no ; then
  AC_MSG_WARN([plumed executable will not run on this machine])
  AC_MSG_WARN([to patch an MD code use 'plumed-patch'])
elif test $mpi_found == ok ; then
  if test $program_can_run_mpi == no ; then
    AC_MSG_WARN([plumed executable will not run on this machine])
    AC_MSG_WARN([unless you invoke it as 'plumed --no-mpi'])
    AC_MSG_WARN([all command line tools are thus available as 'plumed --no-mpi name-of-the-tool'])
    AC_MSG_WARN([e.g. 'plumed --no-mpi driver'])
    AC_MSG_WARN([to patch an MD code use 'plumed --no-mpi patch'])
    AC_MSG_WARN([(notice that MPI will be available anyway in the patched code)])
  fi
fi

AC_SUBST(build_dir)
build_dir=`pwd`



# This is to replace tags in Makefile.conf.in
# saving the result to Makefile.conf
AC_CONFIG_FILES([Makefile.conf sourceme.sh])
# This is to have the stamp-h file up to date
# The date of this file keeps track of when Makefile.conf and sourceme.sh have been updated
AC_CONFIG_FILES([stamp-h], [echo timestamp > stamp-h])
AC_OUTPUT

