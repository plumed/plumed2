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
    
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-$csys-x86_64.sh -O /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -f -p $CONDA_HOME
export PATH="$CONDA_HOME/bin:$PATH"
conda config --set always_yes yes --set changeps1 no
conda config --set anaconda_upload no # not automatically at least
conda update -q conda
conda info -a
conda install conda-build conda-verify anaconda-client

conda-build recipe

ls -l $CONDA_HOME/conda-bld/
ls -l $CONDA_HOME/conda-bld/$TRAVIS_OS_NAME-64


# And now upload if desired
# conda upload $CONDA_HOME/conda-bld/linux-64/*.tar.bz2

# https://gist.github.com/zshaheen/fe76d1507839ed6fbfbccef6b9c13ed9

# https://conda.io/docs/user-guide/tasks/use-conda-with-travis-ci.html


# One could play with this, but perhaps best not to spam the channel repository
export VERSION=`date +%Y.%m.%d`

# Disabled because it fails
if [[ -n "$CONDA_UPLOAD_TOKEN" ]]; then
    USER=plumed  # the conda channel
    anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l testing \
	     $CONDA_HOME/conda-bld/$TRAVIS_OS_NAME-64/plumed*.tar.bz2 --force
fi



