#!/bin/bash

# Anywhere but outside of the repository
export CONDA_HOME=/var/tmp/miniconda

wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -f -p $CONDA_HOME
export PATH="$CONDA_HOME/bin:$PATH"
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
conda install conda-build

conda-build recipe

ls -lR $CONDA_HOME/conda-bld

# And now upload if desired
# conda upload $CONDA_HOME/conda-bld/linux-64/*.tar.bz2


# https://conda.io/docs/user-guide/tasks/use-conda-with-travis-ci.html

