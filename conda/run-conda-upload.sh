#!/bin/bash

# Anywhere but outside of the repository
export CONDA_HOME=/var/tmp/miniconda

if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    csys=Linux
elif [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    csys=MacOSX
else
    echo "Unsupported system $TRAVIS_OS_NAME"
    exit 1
fi
    
export PATH="$CONDA_HOME/bin:$PATH"

# Disabled because it fails
if [[ -n "$CONDA_UPLOAD_TOKEN" ]]; then
    USER=plumed  # the conda channel
    anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l testing \
	     $CONDA_HOME/conda-bld/$TRAVIS_OS_NAME-64/plumed*.tar.bz2 --force
fi



