# The recipe supplies compiler('cuda') in build, with the selected CUDA variant.
CONDA_HOST_REQUIREMENTS = cuda-cudart-dev cuda-cccl  # [linux]
CONDA_RUN_REQUIREMENTS = cuda-cudart  # [linux]
