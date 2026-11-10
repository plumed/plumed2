# Building this plugin requires the NVIDIA HPC SDK compiler.
# It recompiles shared src code with OpenACC; the base package does not need it.
CONDA_UNAVAILABLE_LIBRARIES = NVHPC
