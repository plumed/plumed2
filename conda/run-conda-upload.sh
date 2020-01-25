#!/bin/bash

# Anywhere but outside of the repository
export CONDA_HOME=/var/tmp/miniconda

export PATH="$CONDA_HOME/bin:$PATH"

CONDA_USER=plumed  # the conda channel
if test -n "$TRAVIS_REPO_SLUG" ; then
  CONDA_USER="${TRAVIS_REPO_SLUG%/*}"
fi

if [[ -n "$CONDA_UPLOAD_TOKEN" ]]; then
    anaconda -t $CONDA_UPLOAD_TOKEN upload -u $CONDA_USER -l $CONDA_LABEL $CONDA_HOME/conda-bld/$TRAVIS_OS_NAME-64/plumed*.tar.bz2 --force
    anaconda -t $CONDA_UPLOAD_TOKEN upload -u $CONDA_USER -l $CONDA_LABEL $CONDA_HOME/conda-bld/$TRAVIS_OS_NAME-64/py-plumed*.tar.bz2 --force
fi



