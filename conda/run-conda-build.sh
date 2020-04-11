#!/bin/bash

set -e 

export CPU_COUNT=4
conda-build -c conda-forge plumed
conda-build -c conda-forge py-plumed

ls -l $CONDA_PREFIX/conda-bld/
ls -l $CONDA_PREFIX/conda-bld/$TRAVIS_OS_NAME-64

