#! /usr/bin/env bash

set -e

export CPU_COUNT=4
echo "#Building plumed for conda"
conda-build -c conda-forge plumed
echo "#Building py-plumed for conda"
conda-build -c conda-forge py-plumed

ls -l $CONDA_PREFIX/conda-bld/
ls -l $CONDA_PREFIX/conda-bld/*
