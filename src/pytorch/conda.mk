# CPU packaging first; CUDA variants are a separate packaging dimension.
CONDA_HOST_REQUIREMENTS = "libtorch * cpu_generic*"
CONDA_RUN_REQUIREMENTS = "libtorch * cpu_generic*"
