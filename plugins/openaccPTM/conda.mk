# NVHPC is an external build tool, never a dependency of the main kernel.
# The recipe statically links its runtime and includes the vendor notices.
CONDA_UNAVAILABLE_LIBRARIES = NVHPC
CONDA_HOST_REQUIREMENTS = liblapack
CONDA_RUN_REQUIREMENTS = liblapack
