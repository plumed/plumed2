# Optional native libraries used by this built-in module.
# The base recipe independently selects the features enabled in the kernel.
CONDA_HOST_REQUIREMENTS = "hdf5 * nompi_*" zlib liblapack
CONDA_RUN_REQUIREMENTS = "hdf5 * nompi_*" zlib liblapack
CONDA_UNAVAILABLE_LIBRARIES = NVHPC
