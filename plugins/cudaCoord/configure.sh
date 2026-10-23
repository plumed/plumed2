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
  grep -eDYNAMIC_LIBS -eLDFLAGS -eCXXFLAGS Make.tmp -v
} >Makefile.conf

rm Make.tmp
#tested with nvcc with :"Build cuda_11.7.r11.7/compiler.31442593_0"
#and tested with nvcc with :"Build cuda 11.8
