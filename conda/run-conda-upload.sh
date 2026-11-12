#! /usr/bin/env bash

CONDA_USER=plumed  # the conda channel

if [[ -n "$CONDA_UPLOAD_TOKEN" ]]; then
    anaconda -t $CONDA_UPLOAD_TOKEN upload -u $CONDA_USER -l $CONDA_LABEL $CONDA_PREFIX/conda-bld/*/plumed*.conda --force
    anaconda -t $CONDA_UPLOAD_TOKEN upload -u $CONDA_USER -l $CONDA_LABEL $CONDA_PREFIX/conda-bld/*/py-plumed*.conda --force
fi



