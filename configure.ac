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
# Here we define a few useful macros

# PLUMED_CONFIG_ENABLE(variablename,optionname,doc,default)
# notice that variablename and optionname are likely identical,
# they just need to be different with optionname contains a "-"
# (not allowed in shell variable names)
AC_DEFUN([PLUMED_CONFIG_ENABLE], [
$1=
AC_ARG_ENABLE([$2],
  AS_HELP_STRING([--enable-$2], [enable $3, default: $4]),
  [case "${enableval}" in
             (yes) $1=true ;;
             (no)  $1=false ;;
             (*)   AC_MSG_ERROR([wrong argument to --enable-$2]) ;;
  esac],
  [case "$4" in
             (yes) $1=true ;;
             (no)  $1=false ;;
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

# PLUMED_CHECK_PACKAGE(header,function,define[,library])
# first we check if the header is present. if so, we check if the given function can be found.
# if it cannot be found and the fourth argument (library) has been supplied, we look for it
# in the library. finally, we set the define flag
AC_DEFUN([PLUMED_CHECK_PACKAGE], [
    found=ko
    AC_CHECK_HEADER( [$1],[
      AC_CHECK_FUNC(   [$2], [found=ok],
        m4_ifnblank([$4],[AC_CHECK_LIB(    [$4], [$2], [LIBS="-l$4 $LIBS" && found=ok]) ])
      )
    ])
    if test $found == ok ; then
      AC_DEFINE([$3])
    else
       AC_MSG_WARN([cannot enable $3]) 
    fi
])

##################################################################

# set enable flags for ./configure
PLUMED_CONFIG_ENABLE([debug],[debug],[debugging],[no])
PLUMED_CONFIG_ENABLE([basic_warnings],[basic-warnings],[basic warnings],[yes])
PLUMED_CONFIG_ENABLE([fussy],[fussy],[fussy warnings],[no])
PLUMED_CONFIG_ENABLE([debug_glibcxx],[debug-glibcxx],[enable boundary check],[no])
PLUMED_CONFIG_ENABLE([shared],[shared],[shared libs],[yes])
PLUMED_CONFIG_ENABLE([cxx_exceptions],[cxx-exceptions],[c++ exceptions],[no])
PLUMED_CONFIG_ENABLE([mpi],[mpi],[search for mpi],[yes])
PLUMED_CONFIG_ENABLE([openmp],[openmp],[search for openmp],[no])
PLUMED_CONFIG_ENABLE([external_lapack],[external-lapack],[search for external lapack],[yes])
PLUMED_CONFIG_ENABLE([molfile_plugins],[molfile-plugins],[use molfile_plugins],[yes])
PLUMED_CONFIG_ENABLE([external_molfile_plugins],[external-molfile-plugins],[search for external molfile_plugins],[yes])
PLUMED_CONFIG_ENABLE([matheval],[matheval],[search for matheval],[yes])
PLUMED_CONFIG_ENABLE([zlib],[zlib],[search for zlib],[yes])
PLUMED_CONFIG_ENABLE([clock_gettime],[clock-gettime],[search for gettime],[yes])
PLUMED_CONFIG_ENABLE([gettimeofday],[gettimeofday],[search for gettimeofday],[yes])
PLUMED_CONFIG_ENABLE([cregex],[cregex],[search for C regular expression],[yes])
PLUMED_CONFIG_ENABLE([dlopen],[dlopen],[search for dlopen],[yes])
PLUMED_CONFIG_ENABLE([execinfo],[execinfo],[search for execinfo],[yes])
PLUMED_CONFIG_ENABLE([almost],[almost],[search for almost],[no])
PLUMED_CONFIG_ENABLE([gsl],[gsl],[search for gsl],[no])
PLUMED_CONFIG_ENABLE([xdrfile],[xdrfile],[search for xdrfile],[yes])


# by default use -O flag
# we override the autoconf default (-g) because in release build we do not want to
# include symbol information (obj files are huge)
if test -z "$CXXFLAGS"
then
  CXXFLAGS=-O
fi

# this is list if copied from a generated ./configure script
# the reason why I put it here explicitly is that I want a
# correct fallback list to be searched after preferred choices
# (e.g. MPI) have been checked.
compilers="g++ cxx icpc gpp aCC CC cxx cc++ cl.exe FCC KCC RCC xlC_r xlC"

# if searching for MPI, try first mpi-like compilers
if test $mpi == true ; then
compilers="mpic++ mpicxx mpiicpc openmpic++ openmpicxx $compilers"
fi

# do the actual search
AC_PROG_CXX([$compilers])

# for c, the standard list is used
AC_PROG_CC

# also setup Fortran compiler
# this is optional, and can be used in the late part of this 
# script to verify that fortran can indeed link properly the
# a c++ library
AC_PROG_FC

# we use C++ for all the autoconf tests
AC_LANG(C++)

# log the initial flags
LD="$CXX"
LDSO="$CXX"
AC_MSG_NOTICE([Initial CXX:      $CXX])
AC_MSG_NOTICE([Initial CXXFLAGS: $CXXFLAGS])
AC_MSG_NOTICE([Initial CPPFLAGS: $CPPFLAGS])
AC_MSG_NOTICE([Initial LDFLAGS:  $LDFLAGS])
AC_MSG_NOTICE([Initial LIBS:     $LDFLAGS])
AC_MSG_NOTICE([Initial LD:       $LD])
AC_MSG_NOTICE([Initial LDSO:     $LDSO])
AC_MSG_NOTICE([Initial SOEXT:    $SOEXT])

# check C++ flags
if test $shared == true
then
  PLUMED_CHECK_CXXFLAG([-fPIC])
fi

if test $basic_warnings == true
then
  PLUMED_CHECK_CXXFLAG([-Wall])
  PLUMED_CHECK_CXXFLAG([-pedantic])
  PLUMED_CHECK_CXXFLAG([-ansi])
fi

if test $debug == true
then
  PLUMED_CHECK_CXXFLAG([-g])
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

#### Compulsory libraries ####
# some of them might be made optional if we find that are not available in some system
AC_MSG_NOTICE([Now we will check compulsory headers and libraries])
AC_CHECK_HEADER([dirent.h],     [ ], [AC_MSG_ERROR([compulsory header not found])] )
AC_CHECK_FUNC(  [readdir],      [ ], [AC_MSG_ERROR([compulsory function not found])] )

save_LIBS="$LIBS"
# Then check for blas. This is a bit complicated because we want to allow
# either the version with underscore or the one without
# If they are not found in the normal search path, we repeat the search with -lblas
blas_found=
lapack_found=

if test "$external_lapack" == true ; then

AC_CHECK_FUNC(  [dgemv], [blas_found=nounderscore], [
AC_CHECK_FUNC(  [dgemv_],[blas_found=underscore],   [
AC_CHECK_LIB(   [blas],[dgemv],  [LIBS="-lblas $LIBS"] [blas_found=nounderscore], [
AC_CHECK_LIB(   [blas],[dgemv_], [LIBS="-lblas $LIBS"] [blas_found=underscore]
) ]) ]) ])

if test -n "$blas_found" ; then
# Then we look for lapack using same underscoring
case "$blas_found" in
(underscore) search_for=dsyevr_ ;;
(nounderscore) search_for=dsyevr ;;
esac
AC_CHECK_FUNC( [$search_for], [lapack_found=$blas_found], [
AC_CHECK_LIB(  [lapack],[$search_for], [LIBS="-llapack $LIBS"] [lapack_found=yes]
)
])
fi

fi 

if test -z "$lapack_found" ; then
AC_MSG_WARN([using internal lapack/blas, could be inefficient])
LIBS="$save_LIBS"

AC_DEFINE([__PLUMED_INTERNAL_LAPACK])
AC_DEFINE([__PLUMED_INTERNAL_BLAS])

else

  if test $blas_found == nounderscore
  then
    AC_DEFINE([F77_NO_UNDERSCORE])
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
AC_DEFINE([__PLUMED_HAS_MOLFILE])

found=ko
if test "$external_molfile_plugins" == true ; then
	AC_CHECK_LIB([molfile_plugin],[molfile_dcdplugin_init],
        	[ AC_CHECK_HEADER([libmolfile_plugin.h], [LIBS="-lmolfile_plugin $LIBS"] [found=ok] ) ])
fi
if test $found == ko ; then
	AC_MSG_WARN([using internal molfile_plugins, which only support dcd/xtc/trr/trj/crd files])
	AC_DEFINE([__PLUMED_INTERNAL_MOLFILE_PLUGINS])
fi

fi



# this is special and is also attached to STATIC_LIBS
# this flag should be used also when linking MD engines to allow plumed
# to be loaded later
AC_CHECK_LIB([dl],dlopen, [STATIC_LIBS="-ldl $STATIC_LIBS"] [LIBS="-ldl $LIBS"])

mpi_found=ko
# optional libraries follow
if test $mpi == true ; then
  mpi_found=ok
  PLUMED_CHECK_PACKAGE([mpi.h],[MPI_Init],[__PLUMED_MPI])
else 
  mpi_found=ko
fi
if test $openmp == true ; then
  AC_OPENMP
fi
if test $matheval == true ; then
  PLUMED_CHECK_PACKAGE([matheval.h],[evaluator_create],[__PLUMED_HAS_MATHEVAL],[matheval])
fi
if test $clock_gettime == true ; then
  PLUMED_CHECK_PACKAGE([time.h],[clock_gettime],[__PLUMED_HAS_CLOCK_GETTIME],[rt])
fi
if test $gettimeofday == true ; then
  PLUMED_CHECK_PACKAGE([sys/time.h], [gettimeofday], [__PLUMED_HAS_GETTIMEOFDAY])
fi
if test $cregex == true ; then
  PLUMED_CHECK_PACKAGE([regex.h],[regcomp],[__PLUMED_HAS_CREGEX])
fi
if test $dlopen == true ; then
  PLUMED_CHECK_PACKAGE([dlfcn.h],[dlopen],[__PLUMED_HAS_DLOPEN])
fi
if test $execinfo == true ; then
  PLUMED_CHECK_PACKAGE([execinfo.h],[backtrace],[__PLUMED_HAS_EXECINFO])
fi
if test $zlib == true ; then
  PLUMED_CHECK_PACKAGE([zlib.h],[gzopen],[__PLUMED_HAS_ZLIB],[z])
fi
if test $almost == true ; then
  bz2=ko
  AC_CHECK_HEADER([bzlib.h], AC_SEARCH_LIBS([BZ2_bzReadOpen], [bz2], [bz2=ok] ))
  found=ko
  if test $bz2 == ok; then
    AC_CHECK_HEADER([almost.h], AC_SEARCH_LIBS([init_camshift], [Alm], [found=ok] , AC_MSG_NOTICE(["Almost not found"]), [-lsqlite3 -lz -lbz2 -lnbimpl -lshx]))
  else 
    AC_CHECK_HEADER([almost.h], AC_SEARCH_LIBS([init_camshift], [Alm], [found=ok] , AC_MSG_NOTICE(["Almost not found"]), [-lsqlite3 -lz -lnbimpl -lshx]))
  fi
  if test $found == ok; then
    AC_DEFINE([__PLUMED_HAS_ALMOST]) 
    LIBS="$LIBS -lsqlite3 -lz -lbz2 -lnbimpl -lshx"
  else 
    AC_MSG_ERROR([--enable-almost: cannot find almost])
  fi
fi
if test $gsl == true ; then
  found=ko
  AC_CHECK_LIB([gsl],[main],[found=ok],[found=ko],[-lgslcblas])
  if test $found == ok ; then
    AC_DEFINE([__PLUMED_HAS_GSL])
    LIBS="-lgsl -lgslcblas $LIBS"
  else 
    AC_MSG_ERROR([--enable-gsl: cannot find gsl])
  fi
fi
if test $xdrfile == true ; then
  PLUMED_CHECK_PACKAGE([xdrfile/xdrfile_xtc.h],[write_xtc],[__PLUMED_HAS_XDRFILE],[xdrfile])
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

if test "$cxx_exceptions" == true ; then
  AC_MSG_NOTICE([Enabling c++ exceptions -D__PLUMED_HAS_EXCEPTIONS])
  AC_DEFINE([__PLUMED_HAS_EXCEPTIONS])
fi

# this is necessary in many MPI implementations
# I leave it by default, since it seems harmless
AC_DEFINE([_REENTRANT])

#### Options for dynamic library to work properly ####
AC_SUBST(SOEXT)
AC_SUBST(LD)
AC_SUBST(LDSO)
# these are libraries that should be linked also to MD engines
AC_SUBST(STATIC_LIBS)

if test "$shared" == true ; then
  case `(uname)` in
  (Darwin)
    AC_MSG_NOTICE([*** Special settings for dynamic libraries on OSX ***])
    AC_MSG_NOTICE([Dynamic library extension is 'dylib'])
    AC_MSG_NOTICE([LDSO needs special flags])
    SOEXT=dylib
    LDSO="$LDSO -undefined suppress -flat_namespace -dynamiclib"
  ;;
  (Linux)
    AC_MSG_NOTICE([*** Special settings for dynamic libraries on Linux ***])
    AC_MSG_NOTICE([Dynamic library extension is 'so'])
    AC_MSG_NOTICE([LDSO and LDFLAGS need special flags])
    SOEXT=so
    LDSO="$LDSO -shared"
    PLUMED_CHECK_LDFLAGS([-rdynamic])
  ;;
  (*)
    AC_MSG_NOTICE([*** Dynamic library only enabled on OSX and Linux ***])
  esac
fi

# check linking of runtime library
if test -n "$SOEXT"
then
  AC_MSG_NOTICE([Using LDSO='$LDSO'])
  AC_MSG_NOTICE([Using LDFLAGS='$LDFLAGS'])
  AC_MSG_CHECKING([whether LDSO can create dynamic libraries])
  rm -f conftest.*
  echo "void f(void){ return;}" > conftest.cpp
  $CXX $CXXFLAGS $CPPFLAGS -c conftest.cpp 1>/dev/null 2>/dev/null
  $LDSO $LDFLAGS conftest.o -o conftest.$SOEXT 1>/dev/null 2>/dev/null
  if test -f conftest.$SOEXT
  then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    AC_MSG_WARN([dynamic library will be disabled])
    SOEXT=
  fi
  rm -f conftest.*
fi
#### Options for dynamic library to work properly ####

### Look for doxygen
AC_CHECK_PROG([doxygen],[doxygen],[found])
if test "$doxygen" == found
then
  doxygen_version=`doxygen --version | awk 'BEGIN{FS="."}{if($1>1 || ($1==1 && $2>=8)) print "ok"}'`
  if test "$doxygen_version" == ok
  then
    AC_MSG_NOTICE([Doxygen version is fine])
  else
    AC_MSG_WARN([Doxygen version is <1.8. You might have problems in generating manuals])
  fi
  AC_CHECK_PROG([dot],[dot],[found])
  if test "$dot" != found
  then
    AC_MSG_WARN([You will not be able to see diagrams in the manual])
  fi
fi

### Look for xxd
AC_CHECK_PROG([xxd],[xxd],[found])
if test "$xxd" != found
then
  AC_MSG_ERROR([xxd should be installed for PLUMED to compile properly])
fi

#### This further tests are required to allow linking with non c++ compiler
AC_MSG_NOTICE([PLUMED seems to be configured properly!])
AC_MSG_NOTICE([**************************])
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
    for testlib in "" -lstdc++ -lc++ ; do
      comment=
      test -n "$testlib" && comment=" with library $testlib"
      AC_MSG_CHECKING([whether $name can link a C++ object$comment])
      $compexe $LDFLAGS $testlib $LIBS conftest.o conftest-main.o -o conftest.exe 1>/dev/null 2>/dev/null
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
      AC_MSG_WARN([Please add c++ library to LIBS - e.g. LIBS=-lstdc++])
    fi
  else
    AC_MSG_NOTICE([$compiler compiler not configured])
  fi
  rm -f conftest.* conftest-main.*
done

if test "$prefix" == NONE
then
  prefix=/usr/local
fi

AC_MSG_NOTICE([**** PLUMED will be installed using prefix $prefix])
AC_MSG_NOTICE([**** You can change this later setting PLUMED_PREFIX environment variable before installing])
AC_MSG_NOTICE([**** Executable will be named 'plumed'])
AC_MSG_NOTICE([**** To add a suffix to this name, set PLUMED_LIBSUFFIX environment variable before installing])

if test $mpi == true; then
  if test $mpi_found == ok; then
    AC_MSG_NOTICE([**** PLUMED will be compiled using MPI])
  else 
    AC_MSG_WARN([**** PLUMED will NOT be compiled using MPI because MPI have not been found!])
  fi
else
    AC_MSG_NOTICE([**** PLUMED will be compiled without MPI])
fi

AC_SUBST(build_dir)
build_dir=`pwd`

# This is to replace tags in Makefile.conf.in
# saving the result to Makefile.conf
AC_CONFIG_FILES([Makefile.conf sourceme.sh])
AC_OUTPUT

