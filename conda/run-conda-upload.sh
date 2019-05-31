#!/bin/bash

set -x

# Anywhere but outside of the repository
export CONDA_HOME=/var/tmp/miniconda

export PATH="$CONDA_HOME/bin:$PATH"

if [[ -n "$CONDA_UPLOAD_TOKEN" ]]; then
    USER=plumed  # the conda channel
    anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l $CONDA_LABEL $CONDA_HOME/conda-bld/$TRAVIS_OS_NAME-64/plumed*.tar.bz2 --force
    anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l $CONDA_LABEL $CONDA_HOME/conda-bld/$TRAVIS_OS_NAME-64/py-plumed*.tar.bz2 --force
fi



