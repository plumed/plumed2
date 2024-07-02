### User Configurable Options

## To set up, follow the steps [1], [2], [3] below

#########
## [1]
##
## Set which compiler to use by defining CCCOM:
## GNU GCC compiler
CCCOM=g++ -m64 -std=c++17 -fconcepts -fPIC

## Clang compiler (good to use on macOS)
## Note: try this if you have trouble compiling 
##       with g++ on macOS Mojave 10.14
#CCCOM=clang++ -std=c++17 -fPIC -Wno-gcc-compat

## GCC On Windows cygwin
## Extra flags are workaround for limitations
## of this compiler with Windows binary format
#CCCOM=g++ -std=c++17 -Wa,-mbig-obj -O2 -fPIC

## Intel C++ Compiler not recommended,
## As of June 2016 it does not fully support C++17

#########

#########
## [2]
## 
## BLAS/LAPACK Related Options
## 
## o The variable PLATFORM can be set to macos,lapack,mkl, or acml.
##   This tells ITensor which style of BLAS/Lapack library to expect,
##   and turns various lines of code on or off inside the files
##   itensor/tensor/lapack_wrap.h and lapack_wrap.cc.
## 
## o BLAS_LAPACK_LIBFLAGS specifies blas or lapack related 
##   flags used during the linking step. For example, 
##   flags of the type:
##   -L/path/to/lapack/lib -llapack -lblas
##   though the names of the static libraries (the files referred
##   to by the -l flags) can be highly variable - see examples below.
##
## o BLAS_LAPACK_INCLUDEFLAGS are blas or lapack related flags 
##   needed during compilation. It may include flags such as 
##   -I/path/to/lapack/include
##   where "include" is a folder containing .h header files.
##

##
## Mac OSX system
##

#PLATFORM=macos
#BLAS_LAPACK_LIBFLAGS=-framework Accelerate

##
## Example using a C interface to LAPACK on GNU/LINUX systems
## (Path to lib/ folder may differ on your system)
##

#PLATFORM=lapack
#BLAS_LAPACK_LIBFLAGS=-lpthread -L/usr/lib -lblas -llapack 

##
## Example using the Intel MKL library
## (Path to lib/intel64/ and include/ folders may differ on your system)
##

PLATFORM=mkl
## Example of using a sequential version of MKL. You may want to link to 
## the sequential version of MKL if you are using ITensor's native multithreading, 
## see section [4] for more details.
BLAS_LAPACK_LIBFLAGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
## Use the following line if the library is not in the LD_LIBRARY_PATH
#BLAS_LAPACK_INCLUDEFLAGS=-I/opt/intel/mkl/include

##
## Example using the OpenBLAS library (http://www.openblas.net/)
## (Path to lib/ and include/ folders may differ on your system)
##

#PLATFORM=openblas
#BLAS_LAPACK_LIBFLAGS=-lpthread -L/usr/local/opt/openblas/lib -lopenblas
#BLAS_LAPACK_INCLUDEFLAGS=-I/usr/local/opt/openblas/include -fpermissive -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE

##
## Example using the AMD ACML library
## (Path to lib/ folder may differ on your system)
##

#PLATFORM=acml
#BLAS_LAPACK_LIBFLAGS=-L/opt/acml5.1.0/gfortran64/lib -lacml -lgfortran -lpthread


#########
## [3]
##
## If you want ITensor to be compiled with HDF5 support
## uncomment the line below, and edit it to be the
## location in which the HDF5 libraries and headers are
## installed. 
##
## To determine the prefix on your system, use the command:
## h5cc -show
##
## For more help, visit: 
## http://itensor.org/docs.cgi?vers=cppv3&page=install
##
HDF5_PREFIX=/software/hdf5-1.12.0-el8-x86_64


#########
## [4]
##
## ITensor provides some native OpenMP multithreading.
## Currently, block sparse tensor contractions support
## optional multithreading, so any ITensor code using QN
## conservation may benefit.
##
## If you want to enable OpenMP multithreading in ITensor,
## uncomment the line below. You must have OpenMP installed 
## on your system.
##
## You should set the environment variable OMP_NUM_THREADS
## before running your code to set the number of threads,
## for example using the command:
## export OMP_NUM_THREADS=8
##
## We also recommend setting the number of BLAS/LAPACK threads
## to 1, since BLAS/LAPACK multithreading may compete with ITensor's
## native multithreading.
##
## If you are compiling ITensor with Intel MKL you can set the
## environment variable:
##
## export MKL_NUM_THREADS=1
##
## or directly link to MKL's sequential library (see section [2] above).
##
## If you are compiling ITensor with OpenBLAS, you can set the
## environment variable:
##
## export OPENBLAS_NUM_THREADS=1
##
## For more help, visit: 
## http://itensor.org/docs.cgi?vers=cppv3&page=install
##
ITENSOR_USE_OMP=1


#########
## [5]
##
## This step is optional, but if you wish to customize the flags
## used to compile optimized and debug code, you can do so here.

## Flags to give the compiler for "release mode"
OPTIMIZATIONS=-O2 -DNDEBUG -Wall -Wno-unknown-pragmas

## Flags to give the compiler for "debug mode"
DEBUGFLAGS=-DDEBUG -g -Wall -Wno-unknown-pragmas -pedantic
#
## Set this to 1 if you want ITensor to also build dynamic libraries
## These can be faster to link and give smaller binary sizes
## You may need to set you LD_LIBRARY_PATH to include the ITensor lib/
## folder in order to link with the dynamic libraries
ITENSOR_MAKE_DYLIB=0


###
### Other Makefile variables defined for convenience.
### Not necessary to modify these for most cases.
###

PREFIX=$(THIS_DIR)

ITENSOR_LIBDIR=$(PREFIX)/lib
ITENSOR_INCLUDEDIR=$(PREFIX)

ITENSOR_LIBNAMES=itensor
ITENSOR_LIBFLAGS=$(patsubst %,-l%, $(ITENSOR_LIBNAMES))
ITENSOR_LIBFLAGS+= $(BLAS_LAPACK_LIBFLAGS)
ITENSOR_LIBGFLAGS=$(patsubst %,-l%-g, $(ITENSOR_LIBNAMES))
ITENSOR_LIBGFLAGS+= $(BLAS_LAPACK_LIBFLAGS)
ITENSOR_LIBS=$(patsubst %,$(ITENSOR_LIBDIR)/lib%.a, $(ITENSOR_LIBNAMES))
ITENSOR_GLIBS=$(patsubst %,$(ITENSOR_LIBDIR)/lib%-g.a, $(ITENSOR_LIBNAMES))

ITENSOR_INCLUDEFLAGS=-I'$(ITENSOR_INCLUDEDIR)' $(BLAS_LAPACK_INCLUDEFLAGS)

ifdef HDF5_PREFIX
ITENSOR_USE_HDF5 = 1
ITENSOR_INCLUDEFLAGS += -I$(HDF5_PREFIX)/include -DITENSOR_USE_HDF5
ITENSOR_LIBFLAGS += -L$(HDF5_PREFIX)/lib -lhdf5 -lhdf5_hl
ITENSOR_LIBGFLAGS += -L$(HDF5_PREFIX)/lib -lhdf5 -lhdf5_hl
endif

ifndef CCCOM
$(error Makefile variable CCCOM not defined in options.mk; please define it.)
endif

ifdef ITENSOR_USE_OMP
ITENSOR_INCLUDEFLAGS += -DITENSOR_USE_OMP -fopenmp
ITENSOR_LIBFLAGS += -fopenmp
ITENSOR_LIBGFLAGS += -fopenmp
endif

CCFLAGS=-I. $(ITENSOR_INCLUDEFLAGS) $(OPTIMIZATIONS) -Wno-unused-variable
CCGFLAGS=-I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L'$(ITENSOR_LIBDIR)' $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L'$(ITENSOR_LIBDIR)' $(ITENSOR_LIBGFLAGS)

## Determine shared library extension
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  DYLIB_EXT ?= dylib
  DYLIB_FLAGS ?= -dynamiclib
else
  DYLIB_EXT ?= so
  DYLIB_FLAGS ?= -shared
endif
