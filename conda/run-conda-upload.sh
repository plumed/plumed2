#! /usr/bin/env bash
set -euo pipefail

: "${CONDA_UPLOAD_TOKEN:?CONDA_UPLOAD_TOKEN is required for publication}"
: "${CONDA_LABEL:?CONDA_LABEL is required for publication}"
CONDA_USER=${CONDA_USER:-plumed}
if [[ -n ${TRAVIS_REPO_SLUG:-} ]]; then
  CONDA_USER=${TRAVIS_REPO_SLUG%/*}
fi
build_root=${CONDA_BLD_PATH:-${CONDA_PREFIX:?}/conda-bld}
version=${VERSION:-*}
shopt -s nullglob
packages=("$build_root"/*/plumed*-"$version"-*.conda "$build_root"/*/py-plumed-"$version"-*.conda)
# An unset VERSION keeps the historical all-version upload behavior.
if [[ $version == '*' ]]; then
  packages=("$build_root"/*/plumed*.conda "$build_root"/*/py-plumed*.conda)
fi
if (( ${#packages[@]} == 0 )); then
  echo "No PLUMED package archives found for version $version in $build_root" >&2
  exit 1
fi
# Read by anaconda-client; keep credentials out of command arguments.
export ANACONDA_API_TOKEN=$CONDA_UPLOAD_TOKEN
arguments=(upload --user "$CONDA_USER" --label "$CONDA_LABEL")
# Dated snapshots have unique build numbers and must not replace an earlier set.
if [[ $CONDA_LABEL != test ]]; then
  arguments+=(--force)
fi
anaconda "${arguments[@]}" "${packages[@]}"
