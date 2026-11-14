# Options: see conda/module-options.md (repository root).
# NVHPC is an external build tool, never a dependency of the main kernel.
# The recipe statically links its runtime and includes the vendor notices.
CONDA_UNAVAILABLE_LIBRARIES = NVHPC
CONDA_HOST_REQUIREMENTS = libopenblas
CONDA_RUN_REQUIREMENTS = libopenblas
# OpenBLAS supplies the directly imported LAPACK routines and its own SONAME.

CONDA_PLATFORMS = linux-64
CONDA_VARIANT_nvhpc_version = 25.11
CONDA_VARIANT_openacc_cuda_version = 13.0
CONDA_ZIP_KEYS = nvhpc_version,openacc_cuda_version
CONDA_VARIANT_c_compiler_version = 14
CONDA_VARIANT_cxx_compiler_version = 14
CONDA_VARIANT_c_stdlib = sysroot
CONDA_VARIANT_c_stdlib_version = 2.28
CONDA_TEST_REQUIREMENTS = python
CONDA_TEST_SOURCE_FILES = extensions/tests/openacc.py plugins/openaccPTM/regtest/targetGPU/rt-GPU-DISTANCE plugins/openaccPTM/regtest/targetSAXS/rt-GPU-saxs-chi2-deriv regtest/isdb/rt-saxs-chi2-deriv plugins/openaccPTM/regtest/trajectories/isomers.xyz regtest/scripts
CONDA_TEST_COMMANDS = 'python extensions/tests/openacc.py --prefix "$PREFIX" --device host'
