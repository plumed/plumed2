compile := nvcc -g -ccbin mpic++ -c -I"$PLUMED_INCLUDEDIR/$PLUMED_PROGRAM_NAME"/ -DPACKAGE_NAME="PLUMED" -DPACKAGE_TARNAME="plumed" -DPACKAGE_VERSION="2" -DPACKAGE_STRING="PLUMED 2" -DPACKAGE_BUGREPORT="" -DPACKAGE_URL="" -D__PLUMED_LIBCXX11=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -D__PLUMED_HAS_EXTERNAL_BLAS=1 -D__PLUMED_HAS_EXTERNAL_LAPACK=1 -D__PLUMED_HAS_MOLFILE_PLUGINS=1 -D__PLUMED_HAS_MPI=1 -D__PLUMED_HAS_ASMJIT=1 -D__PLUMED_HAS_CREGEX=1 -D__PLUMED_HAS_DLOPEN=1 -D__PLUMED_HAS_RTLD_DEFAULT=1 -D__PLUMED_HAS_CHDIR=1 -D__PLUMED_HAS_SUBPROCESS=1 -D__PLUMED_HAS_GETCWD=1 -D__PLUMED_HAS_POPEN=1 -D__PLUMED_HAS_EXECINFO=1 -D__PLUMED_HAS_ZLIB=1 -D__PLUMED_HAS_GSL=1 -D__PLUMED_HAS_FFTW=1 -DNDEBUG=1 -D_REENTRANT=1 -O3 -Xcompiler -fPIC -Xcompiler -Wall -Xcompiler -std=c++11 -Xcompiler -fno-gnu-unique -Xcompiler -fopenmp
link    := nvcc -shared -lstdc++ -lfftw3 -lgsl -lgslcblas -lz -ldl -llapack -lblas -Xcompiler -rdynamic -Xlinker -Bsymbolic -Xcompiler -fopenmp /u/d/drapetti/scratch/installs/plumedCuda/lib/libplumedKernel.so
all: Coordination.so

%.o: %.cu ndReduction.h cudaHelpers.cuh
	$(compile) -o $@  $<

Coordination.so: ndReduction.o Coordination.o 
	$(link) -o $@ $^

clean:
	rm Coordination.so ndReduction.o 
