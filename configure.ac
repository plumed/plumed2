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
PLUMED_CONFIG_ENABLE([doc],[doc],[documentation],[yes])
PLUMED_CONFIG_ENABLE([pdfdoc],[pdfdoc],[pdf version of the manual],[no])
PLUMED_CONFIG_ENABLE([debug],[debug],[debugging],[no])
PLUMED_CONFIG_ENABLE([gcov],[gcov],[gcov to estimate code coverage],[no])
PLUMED_CONFIG_ENABLE([basic_warnings],[basic-warnings],[basic warnings],[yes])
PLUMED_CONFIG_ENABLE([fussy],[fussy],[fussy warnings],[no])
PLUMED_CONFIG_ENABLE([debug_glibcxx],[debug-glibcxx],[enable boundary check],[no])
PLUMED_CONFIG_ENABLE([shared],[shared],[shared libs],[yes])
PLUMED_CONFIG_ENABLE([dependency_tracking],[dependency-tracking],[dependency tracking],[yes])
PLUMED_CONFIG_ENABLE([cxx_exceptions],[cxx-exceptions],[c++ exceptions],[no])
PLUMED_CONFIG_ENABLE([ld_r],[ld-r],[group object files],[yes])
PLUMED_CONFIG_ENABLE([mpi],[mpi],[search for mpi],[yes])
PLUMED_CONFIG_ENABLE([external_lapack],[external-lapack],[search for external lapack],[yes])
PLUMED_CONFIG_ENABLE([external_blas],[external-blas],[search for external blas],[yes])
PLUMED_CONFIG_ENABLE([molfile_plugins],[molfile-plugins],[use molfile_plugins],[yes])
PLUMED_CONFIG_ENABLE([external_molfile_plugins],[external-molfile-plugins],[search for external molfile_plugins],[yes])
PLUMED_CONFIG_ENABLE([matheval],[matheval],[search for matheval],[yes])
PLUMED_CONFIG_ENABLE([zlib],[zlib],[search for zlib],[yes])
PLUMED_CONFIG_ENABLE([clock_gettime],[clock-gettime],[search for gettime],[yes])
PLUMED_CONFIG_ENABLE([gettimeofday],[gettimeofday],[search for gettimeofday],[yes])
PLUMED_CONFIG_ENABLE([readdir_r],[readdir-r],[search for readdir_r (threadsafe)],[yes])
PLUMED_CONFIG_ENABLE([cregex],[cregex],[search for C regular expression],[yes])
PLUMED_CONFIG_ENABLE([dlopen],[dlopen],[search for dlopen],[yes])
PLUMED_CONFIG_ENABLE([execinfo],[execinfo],[search for execinfo],[yes])
PLUMED_CONFIG_ENABLE([gsl],[gsl],[search for gsl],[yes])
PLUMED_CONFIG_ENABLE([xdrfile],[xdrfile],[search for xdrfile],[yes])
PLUMED_CONFIG_ENABLE([boost_graph],[boost_graph],[search for boost graph],[no])
PLUMED_CONFIG_ENABLE([fftw],[fftw],[search for fftw],[no])

AC_ARG_VAR(SOEXT,[extension of dynamic libraries (so/dylib)])
AC_ARG_VAR(STATIC_LIBS,[variables that should be linked statically directly to MD code - configure will add here -ldl if necessary ])
AC_ARG_VAR(LDSO,[command for linking shared library - configure will use CXX plus the proper flags ])


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
AC_MSG_NOTICE([Initial CXX:         $CXX])
AC_MSG_NOTICE([Initial CXXFLAGS:    $CXXFLAGS])
AC_MSG_NOTICE([Initial CPPFLAGS:    $CPPFLAGS])
AC_MSG_NOTICE([Initial CFLAGS:      $CFLAGS])
AC_MSG_NOTICE([Initial LDFLAGS:     $LDFLAGS])
AC_MSG_NOTICE([Initial LIBS:        $LIBS])
AC_MSG_NOTICE([Initial STATIC_LIBS: $STATIC_LIBS])
AC_MSG_NOTICE([Initial LD:          $LD])
AC_MSG_NOTICE([Initial LDSO:        $LDSO])
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
  PLUMED_CHECK_CXXFLAG([-ansi])
fi

if test $debug == true
then
  PLUMED_CHECK_CXXFLAG([-g])
fi

if test $gcov == true
then
  PLUMED_CHECK_CXXFLAG([--coverage])
  CFLAGS="$CFLAGS --coverage"
  LDFLAGS="$CFLAGS --coverage"
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

AC_MSG_CHECKING([whether $CXX supports explicit])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([ class A{explicit A(){}};])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no]) ;
   AC_DEFINE([explicit],[])])

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

save_LIBS="$LIBS"
# Then check for blas. This is a bit complicated because we want to allow
# either the version with underscore or the one without
# If they are not found in the normal search path, we repeat the search with -lblas
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
AC_CHECK_FUNC(  [dgemv], [blas_found=nounderscore], [
AC_CHECK_FUNC(  [dgemv_],[blas_found=underscore],   [
AC_CHECK_LIB(   [blas],[dgemv],  [LIBS="-lblas $LIBS"] [blas_found=nounderscore], [
AC_CHECK_LIB(   [blas],[dgemv_], [LIBS="-lblas $LIBS"] [blas_found=underscore]
) ]) ]) ])
fi

# if not found, then use internal lapack and blas
if test -z "$blas_found" ; then
AC_MSG_WARN([using internal lapack and blas, could be inefficient])
LIBS="$save_LIBS"

fi

# if found, also look for external lapack
if test -n "$blas_found" ; then

AC_DEFINE([__PLUMED_HAS_EXTERNAL_BLAS])

save_LIBS="$LIBS"

if test "$external_lapack" == true ; then
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

# if not found, then use internal lapack with external blas
if test -z "$lapack_found" ; then
AC_MSG_WARN([using internal lapack, could be inefficient])
LIBS="$save_LIBS"

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

found=ko
if test "$external_molfile_plugins" == true ; then
	AC_CHECK_LIB([molfile_plugin],[molfile_dcdplugin_init],
        	[ AC_CHECK_HEADER([libmolfile_plugin.h], [LIBS="-lmolfile_plugin $LIBS"] [found=ok] ) ])
fi
if test $found == ko ; then
	AC_MSG_WARN([using internal molfile_plugins, which only support dcd/xtc/trr/trj/crd files])
else
	AC_DEFINE([__PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS])
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
  if test "x$ac_cv_func_MPI_Init" = xyes; then
    mpi_found=ok
  fi
else 
  mpi_found=ko
fi

# search for openmp is automatically disabled by autoconf
# when configuring with --disable-openmp
AC_OPENMP

if test $matheval == true ; then
  PLUMED_CHECK_PACKAGE([matheval.h],[evaluator_create],[__PLUMED_HAS_MATHEVAL],[matheval])
fi
if test $clock_gettime == true ; then
  PLUMED_CHECK_PACKAGE([time.h],[clock_gettime],[__PLUMED_HAS_CLOCK_GETTIME],[rt])
fi
if test $gettimeofday == true ; then
  PLUMED_CHECK_PACKAGE([sys/time.h], [gettimeofday], [__PLUMED_HAS_GETTIMEOFDAY])
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
if test $execinfo == true ; then
  PLUMED_CHECK_PACKAGE([execinfo.h],[backtrace],[__PLUMED_HAS_EXECINFO])
fi
if test $zlib == true ; then
  PLUMED_CHECK_PACKAGE([zlib.h],[gzopen],[__PLUMED_HAS_ZLIB],[z])
fi

save_LIBS="$LIBS"

if test $gsl == true ; then
  found=ko
  AC_CHECK_FUNC(  [cblas_dgemv], [found=ok], [
  AC_CHECK_LIB(   [gslcblas],[cblas_dgemv], [LIBS="-lgslcblas $LIBS"] [found=ok]
  ) ])
  if test $found == ok ; then
    found=ko
    AC_CHECK_HEADER(  [gsl/gsl_vector.h], [
      AC_CHECK_FUNC(  [gsl_vector_alloc], [found=ok], [
      AC_CHECK_LIB(   [gsl],[gsl_vector_alloc], [LIBS="-lgsl $LIBS"] [found=ok]
      ) ])
    ])
  fi
  if test $found == ok ; then
    AC_DEFINE([__PLUMED_HAS_GSL])
  else
     LIBS="$save_LIBS"
     AC_MSG_WARN([cannot enable __PLUMED_HAS_GSL])
  fi
fi

if test $xdrfile == true ; then
  PLUMED_CHECK_PACKAGE([xdrfile/xdrfile_xtc.h],[write_xtc],[__PLUMED_HAS_XDRFILE],[xdrfile])
fi
if test $boost_graph == true ; then
  PLUMED_CHECK_PACKAGE([boost/graph/graph_utility.hpp],[exit],[__PLUMED_HAS_BOOST_GRAPH],[boost_graph])
fi
if test $fftw == true ; then
  PLUMED_CHECK_PACKAGE([fftw3.h],[fftw_execute],[__PLUMED_HAS_FFTW],[fftw3])
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

#### This further tests are required to allow linking with non c++ compiler
AC_MSG_NOTICE([PLUMED seems to be configured properly!])
AC_MSG_NOTICE([**************************])
AC_SUBST(LD_RO)
LD_RO=
if test "$ld_r" == true ; then
  AC_MSG_CHECKING([whether C++ objects can be grouped with ld -r])
  
  LD_RO="ld -r -o"
  
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
  
  $LD_RO conftest-both.o conftest-f.o conftest-g.o 1> /dev/null 2> /dev/null
  
  $CXX $CXXFLAGS -o conftest.exe conftest-main.o conftest-both.o 1> /dev/null 2> /dev/null
  
  if test -f conftest.exe
  then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    LD_RO=""
  fi
fi

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
    for testlib in "" -lstdc++ -lc++ -lmpi_cxx ; do
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
      AC_MSG_WARN([Please add c++ library to LIBS])
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

AC_ARG_PROGRAM


# version modified to work in shell script instead of makefile:
program_transform_name_sh=$(echo "${program_transform_name}" | sed 's:\$\$:$:g')
program_name=$(echo plumed | sed "$program_transform_name_sh")
AC_SUBST(program_name)

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
AC_MSG_NOTICE([**** You can change paths later using options to "make install"])
AC_MSG_NOTICE([**** e.g. with "make install prefix=/path"])

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

