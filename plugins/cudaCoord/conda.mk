# Options: see conda/module-options.md (repository root).
# The recipe supplies compiler('cuda') in build, with the selected CUDA variant.
CONDA_HOST_REQUIREMENTS = cuda-cudart-dev cuda-cccl  # [linux]
CONDA_RUN_REQUIREMENTS = cuda-cudart  # [linux]

CONDA_PLATFORMS = linux-64 linux-aarch64
CONDA_VARIANT_cuda_compiler = cuda-nvcc
CONDA_VARIANT_cuda_compiler_version = 12.9 13.0
CONDA_VARIANT_cxx_compiler_version = 14
CONDA_VARIANT_c_stdlib = sysroot
CONDA_VARIANT_c_stdlib_version = 2.28
CONDA_TEST_REQUIREMENTS = python
CONDA_TEST_SOURCE_FILES = extensions/tests/cuda_coordination.py plugins/cudaCoord/regtest/cudatest/rt-double-nopbc plugins/cudaCoord/regtest/cudatest/rt-float-nopbc plugins/cudaCoord/regtest/trajectories/traj-nopbc.xyz regtest/scripts
CONDA_TEST_COMMANDS = 'python extensions/tests/cuda_coordination.py --prefix "$PREFIX" --cuda-version "{{ cuda_compiler_version }}"'
