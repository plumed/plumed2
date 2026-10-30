#!/bin/bash
#formatted with shfmt -ci -s -bn -i 2
#checked with shellcheck
unset CXX
unset CC
unset CFLAGS
unset CXXFLAGS
unset LDSHARED
export plumed_program_name=plumed
# plumed.c (generated from plumed.pyx and cplumed.pxd, see "make cythonize" in python/Makefile)
# is compiled directly, so cython is not needed and plumed.c is never modified by the build.
# Export plumed_force_cython=yes before calling this script to compile plumed.pyx instead.
export plumed_include_dir=../src/wrapper
while getopts ":v:k:p:" opt; do
  case $opt in
    v)
      export plumed_version=${OPTARG}
      ;;
    k)
      export plumed_default_kernel=${OPTARG}
      ;;
    p)
      python_bin=${OPTARG}
      ;;
    :) echo "Wrong inputs use -v for version, -k for kernel and -p for python exec" ;;
    \?) echo "Wrong inputs use -v for version, -k for kernel and -p for python exec" ;;
  esac
done
echo "plumed_version=${plumed_version}"
echo "plumed_default_kernel=${plumed_default_kernel}"
echo "python_bin=${python_bin}"
${python_bin} -m build --wheel .
#getting the shared object out of the wheel
unzip ./dist/plumed-*.whl -x "*/*"
#some cleaning
rm -rf dist build
