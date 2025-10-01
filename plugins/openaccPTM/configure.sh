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

#pendantic adds a unuseful FOR EACH line with
#"" warning: style of line directive is a GCC extension"
{
    #    grep CXXFLAGS Make.tmp \
    #        | sed -e 's/-f/-Xcompiler -f/g' \
    #            -e 's/-pedantic//g' \
    #            -e 's/-W/-Xcompiler -W/g'
    #    grep -eDYNAMIC_LIBS -eLDFLAGS Make.tmp \
    #        | sed -e 's/-rdynamic/-Xcompiler -rdynamic/g' \
    #            -e 's/-Wl,/-Xlinker /g' \
    #            -e 's/-f/-Xcompiler -f/g'
    #    #prints the rest of the file
    #    grep -eDYNAMIC_LIBS -eLDFLAGS -eCXXFLAGS Make.tmp -v
    sed -e 's/-fno-gnu-unique//g' Make.tmp
} >Makefile.conf

rm Make.tmp
