#!/usr/bin/env bash
set -euo pipefail
export PLUMED_CONDA_BUILD=1
./configure --with-extension-base="$PREFIX" --enable-extension-modules=all --prefix="$PREFIX"
make
make install
