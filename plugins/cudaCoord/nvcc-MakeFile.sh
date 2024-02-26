#! /usr/bin/env bash

if [[ -z $PLUMED_KERNEL ]]; then
  echo "$(basename $0) can work only if \"PLUMED_KERNEL\" is defined"
  echo "either via module load or sourceme.sh"
fi

{
  plumed config makefile_conf
  echo "PLUMED_INCLUDE=-I$(plumed info --include-dir)"
  echo "PLUMED_KERNEL=-L${PLUMED_KERNEL}"
} >Make.tmp

if [[ ${SILENT_CUDA_COMPILATION} ]]; then
  #-w suppress the warnings
  sed -i -e 's/-Wall//g' \
    -e 's/-c/-w -c /g' Make.tmp
fi

#pendantic adds a unuseful FOR EACH line with
#"" warning: style of line directive is a GCC extension"
{
  grep CXXFLAGS Make.tmp |
    sed -e 's/-f/-Xcompiler -f/g' \
      -e 's/-pedantic//g' \
      -e 's/-W/-Xcompiler -W/g'
  grep -eDYNAMIC_LIBS -eLDFLAGS Make.tmp |
    sed -e 's/-rdynamic/-Xcompiler -rdynamic/g' \
      -e 's/-Wl,/-Xlinker /g' \
      -e 's/-f/-Xcompiler -f/g'
  #prints the rest of the file
  grep -eDYNAMIC_LIBS -eLDFLAGS -eCXXFLAGS Make.tmp -v
} >Makefile.conf

rm Make.tmp
#tested with nvcc with :"Build cuda_11.7.r11.7/compiler.31442593_0"
#and tested with nvcc with :"Build cuda 11.8
