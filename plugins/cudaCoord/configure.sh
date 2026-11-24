#! /usr/bin/env bash

if [[ -z $PLUMED_KERNEL ]]; then
  echo "$(basename $0) can work only if \"PLUMED_KERNEL\" is defined"
  echo "either via module load or sourceme.sh"
fi

{
  plumed config makefile_conf
  echo "PLUMED_INCLUDE=-I$(plumed info --include-dir)"
  #link the plugin directly to the kernel (and remember where it is), so that it
  #can be LOADed also when the kernel has been loaded with RTLD_LOCAL
  #(the default for runtime loading since PLUMED 2.10, e.g. GROMACS with PLUMED_KERNEL)
  echo "PLUMED_KERNEL=${PLUMED_KERNEL} -Xlinker -rpath=$(dirname "${PLUMED_KERNEL}")"
} >Make.tmp

if [[ ${SILENT_CUDA_COMPILATION} ]]; then
  #-w suppress the warnings
  sed -i -e 's/-Wall//g' \
    -e 's/-c/-w -c /g' Make.tmp
fi

if [[ ${USE_HIP} ]]; then
  #hipcc (clang) takes the flags of the host compiler as they are, except
  #-fconcepts and -fno-gnu-unique, that it does not know.
  #-fopenmp is not given to the linker: it would add the OpenMP runtime of clang
  #to the one of the kernel, and the plugin does not call OpenMP.
  #hipcc does not use the host compiler, that may be an MPI wrapper, so the
  #directory of mpi.h is asked to the host compiler
  {
    grep CXXFLAGS Make.tmp \
      | sed -E -e 's/-fconcepts\>//g' \
        -e 's/-fno-gnu-unique\>//g'
    grep -eDYNAMIC_LIBS -eLDFLAGS Make.tmp \
      | sed -E -e 's/ -fopenmp( |$)/\1/g'
    grep -eDYNAMIC_LIBS -eLDFLAGS -eCXXFLAGS Make.tmp -v
    if grep -q -e '-D__PLUMED_HAS_MPI' Make.tmp; then
      read -r -a cxx <<<"$(sed -n 's/^CXX=//p' Make.tmp)"
      mpih=$(echo '#include <mpi.h>' | "${cxx[@]}" -x c++ -E -H -o /dev/null - 2>&1 \
        | sed -n 's|^\. \(.*/mpi\.h\)$|\1|p')
      if [[ -n $mpih ]]; then
        echo "CPPFLAGS+=-I${mpih%/mpi.h}"
      else
        echo "$(basename "$0"): mpi.h not found by ${cxx[0]}" >&2
      fi
    fi
  } >Makefile.conf
  rm Make.tmp
  exit 0
fi

#-fconcepts (e.g. needed by ITensor in the main build) makes libstdc++ expose C++20
#concepts code that nvcc cannot parse in C++17 mode, so it is removed
#pendantic adds a unuseful FOR EACH line with
#"" warning: style of line directive is a GCC extension"
{
  grep CXXFLAGS Make.tmp \
    | sed -E -e 's/-fconcepts\>//g' \
      -e 's/-f([a-zA-Z-]*)\>/-Xcompiler -f\1/g' \
      -e 's/-gdwarf/-Xcompiler -gdwarf/g' \
      -e 's/-pedantic//g' \
      -e 's/-W/-Xcompiler -W/g'
  grep -eDYNAMIC_LIBS -eLDFLAGS Make.tmp \
    | sed -e 's/-rdynamic/-Xcompiler -rdynamic/g' \
      -e 's/-Wl,/-Xlinker /g' \
      -e 's/-f/-Xcompiler -f/g'
  #prints the rest of the file
  #CXX is passed to nvcc -ccbin, that takes a compiler and not a launcher:
  #ccache (e.g. PLUMED configured with CXX="ccache mpic++") is removed
  grep -eDYNAMIC_LIBS -eLDFLAGS -eCXXFLAGS Make.tmp -v \
    | sed -E -e 's/^CXX=ccache[[:space:]]+/CXX=/'
} >Makefile.conf

rm Make.tmp
#tested with nvcc with :"Build cuda_11.7.r11.7/compiler.31442593_0"
#and tested with nvcc with :"Build cuda 11.8
