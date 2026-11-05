#!/usr/bin/env bash
set -e

source ../sourceme.sh
../src/lib/plumed --no-mpi gen_vim --outdir . --template-dir templates
