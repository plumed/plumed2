# Build against serial HDF5, but allow ABI-compatible MPI variants at runtime.
# Retain the main package's configured FFTW capability and actual OpenBLAS provider.
CONDA_HOST_REQUIREMENTS = "hdf5 * nompi_*" zlib libopenblas "fftw * nompi_*"
CONDA_RUN_REQUIREMENTS = hdf5 zlib libopenblas "fftw * nompi_*"
