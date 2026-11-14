# Options: see conda/module-options.md (repository root).
# CPU packaging first; CUDA variants are a separate packaging dimension.
CONDA_HOST_REQUIREMENTS = "libtorch * cpu_generic*"
CONDA_RUN_REQUIREMENTS = "libtorch * cpu_generic*"

CONDA_TEST_REGTESTS = pytorch/rt-pytorch_model_strided
