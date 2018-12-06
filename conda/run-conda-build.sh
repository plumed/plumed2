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
conda update -q conda
conda info -a
conda install conda-build conda-verify

conda-build recipe

ls -l $CONDA_HOME/conda-bld/
ls -lR $CONDA_HOME/conda-bld/*-64  $CONDA_HOME/conda-bld/noarch

# And now upload if desired
# conda upload $CONDA_HOME/conda-bld/linux-64/*.tar.bz2


# https://conda.io/docs/user-guide/tasks/use-conda-with-travis-ci.html

