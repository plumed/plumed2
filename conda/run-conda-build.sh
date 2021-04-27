#!/bin/bash

set -e 

if [[ "$(uname)" == Linux ]]; then
    OS_NAME=linux
elif [[ "$(uname)" == "Darwin" ]]; then
    OS_NAME=osx
else
    echo "Unsupported system $(uname)"
    exit 1
fi
    
export CPU_COUNT=4
conda-build -c conda-forge plumed
conda-build -c conda-forge py-plumed

ls -l $CONDA_PREFIX/conda-bld/
ls -l $CONDA_PREFIX/conda-bld/$OS_NAME-64

