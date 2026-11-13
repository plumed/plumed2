# Options: see conda/module-options.md (repository root).
# EMMIVOX uses LibTorch directly; it does not link the PLUMED pytorch module.
CONDA_HOST_REQUIREMENTS = "libtorch * cpu_generic*"
CONDA_RUN_REQUIREMENTS = "libtorch * cpu_generic*"

CONDA_TEST_REGTESTS = isdb-pytorch/rt-emmivox-fixed isdb-pytorch/rt-checkpoint-emmivox
