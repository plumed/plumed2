#!/usr/bin/env bash
# Run already-built packages on NVIDIA hardware; never override driver metadata.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VERSION:?set VERSION to the artifact package version}"
: "${GPU_CONDA_CHANNEL:?set GPU_CONDA_CHANNEL to the indexed artifact directory}"
: "${GPU_TEST_ROOT:?set GPU_TEST_ROOT to a writable test directory}"
unset CONDA_OVERRIDE_CUDA
nvidia-smi
mkdir -p "$GPU_TEST_ROOT"
failed=0
for variant in cuda12.9 cuda13.0 openacc; do
  prefix="$GPU_TEST_ROOT/$variant"
  if [[ $variant == cuda* ]]; then
    build_variant=${variant//./}
    spec="plumed-ext-cuda-coordination=$VERSION=${build_variant}_*"
  else
    spec="plumed-ext-openacc-ptm=$VERSION"
  fi
  if ! conda create -y -p "$prefix" --override-channels -c "$GPU_CONDA_CHANNEL" -c conda-forge "plumed=$VERSION" "$spec" python=3.13 > "$GPU_TEST_ROOT/$variant.log" 2>&1; then
    echo "$variant: installation failed (see $GPU_TEST_ROOT/$variant.log)"
    failed=1
    continue
  fi
  if [[ $variant == cuda* ]]; then
    command=("$prefix/bin/python" extensions/tests/cuda_coordination.py --prefix "$prefix" --cuda-version "${variant#cuda}" --gpu)
  else
    command=("$prefix/bin/python" extensions/tests/openacc.py --prefix "$prefix" --device nvidia)
  fi
  if "${command[@]}" >> "$GPU_TEST_ROOT/$variant.log" 2>&1; then
    echo "$variant: GPU numerical tests passed"
  else
    echo "$variant: GPU numerical tests failed (see $GPU_TEST_ROOT/$variant.log)"
    failed=1
  fi
  cat "$GPU_TEST_ROOT/$variant.log"
done
exit "$failed"
