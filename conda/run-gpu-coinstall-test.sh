#!/usr/bin/env bash
# Native loading/solver validation only: no GPU calculation is executed here.
set -euo pipefail
[[ $(uname) == Linux && $(uname -m) == x86_64 ]] || exit 0
: "${VERSION:?set VERSION to the built package version}"
: "${CONDA_BLD_PATH:?set CONDA_BLD_PATH to the local channel}"
: "${CPU_CONDA_ARTIFACT_DIR:?set CPU_CONDA_ARTIFACT_DIR to the CPU/CUDA artifacts}"
: "${OPENACC_CONDA_ARTIFACT_DIR:?set OPENACC_CONDA_ARTIFACT_DIR to the OpenACC artifacts}"
python - <<'PY'
import os
import shutil
from pathlib import Path

sources = [Path(os.environ[name]) for name in ('CPU_CONDA_ARTIFACT_DIR', 'OPENACC_CONDA_ARTIFACT_DIR')]
version = os.environ['VERSION']
bases = [list(root.rglob(f'plumed-{version}-*.conda')) for root in sources]
assert all(len(files) == 1 for files in bases), bases
assert bases[0][0].read_bytes() == bases[1][0].read_bytes(), 'GPU jobs used different kernel archives'
destination = Path(os.environ['CONDA_BLD_PATH']) / 'linux-64'
destination.mkdir(parents=True, exist_ok=True)
for root in sources:
    for package in root.rglob('*.conda'):
        target = destination / package.name
        if target.exists():
            assert target.read_bytes() == package.read_bytes(), f'Conflicting archive: {package.name}'
        else:
            shutil.copy2(package, target)
print('CPU/CUDA and OpenACC artifacts contain the exact same kernel archive')
PY
conda index "$CONDA_BLD_PATH"
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-gpu-coinstall.XXXXXX")
trap 'rm -rf "$task"' EXIT
prefix="$task/environment"
# Override driver metadata solely to check solver selection and native loading.
export CONDA_OVERRIDE_CUDA=${CONDA_OVERRIDE_CUDA:-$(python conda/gpu-plan.py --max-cuda)} NVCOMPILER_ACC_DEVICE_TYPE=HOST
conda create -y -p "$prefix" --override-channels -c local -c conda-forge "plumed=$VERSION"
cp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
for cuda_version in $(python conda/gpu-plan.py --cuda-versions); do
  conda install -y -p "$prefix" --override-channels -c local -c conda-forge \
    "plumed-ext-cuda-coordination=$VERSION" "cuda-version=$cuda_version" "plumed-ext-openacc-ptm=$VERSION"
  cmp "$prefix/lib/libplumedKernel.so" "$task/kernel-before"
  "$prefix/bin/plumed" modules -q --has cuda_coordination openacc_ptm
  "$prefix/bin/plumed" modules --module cuda_coordination --has-feature cuda -q
  "$prefix/bin/plumed" modules --module openacc_ptm --has-feature openacc -q
done
echo 'GPU packages co-install and load against one unchanged kernel; GPU numerics untested'
