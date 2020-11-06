#!/bin/bash

set -e 

# Anywhere but outside of the repository
export CONDA_HOME=/var/tmp/miniconda

if [[ "$(uname)" == Linux ]]; then
    csys=Linux
    OS_NAME=linux
elif [[ "$(uname)" == "Darwin" ]]; then
    csys=MacOSX
    OS_NAME=osx
    PREV=$(pwd)
    cd /var/tmp
    mkdir MacOSX-SDKs
    cd MacOSX-SDKs
    wget https://github.com/phracker/MacOSX-SDKs/releases/download/10.13/MacOSX10.9.sdk.tar.xz
    tar -xf ./MacOSX10.9.sdk.tar.xz
    rm MacOSX10.9.sdk.tar.xz
    cd $PREV
else
    echo "Unsupported system $(uname)"
    exit 1
fi
    
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-$csys-x86_64.sh -O /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -f -p $CONDA_HOME
export PATH="$CONDA_HOME/bin:$PATH"
conda config --set always_yes yes --set changeps1 no
conda config --set anaconda_upload no # not automatically at least
conda update -q conda
conda info -a
conda install conda-build conda-verify anaconda-client

export CPU_COUNT=4

conda-build -c conda-forge plumed
conda-build -c conda-forge py-plumed

ls -l $CONDA_HOME/conda-bld/
ls -l $CONDA_HOME/conda-bld/$OS_NAME-64

