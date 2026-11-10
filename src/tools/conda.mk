# Build against serial HDF5, but allow ABI-compatible MPI variants at runtime.
# Optional native libraries used by this built-in module.
# The base recipe independently selects the features enabled in the kernel.
CONDA_HOST_REQUIREMENTS = "hdf5 * nompi_*" zlib liblapack
CONDA_RUN_REQUIREMENTS = hdf5 zlib liblapack
# OpenACC helpers are compiled with NVIDIA tools only in plugins/openaccPTM.
