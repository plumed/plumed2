# Build against serial HDF5, but allow ABI-compatible MPI variants at runtime.
# FFTW and the OpenBLAS provider supply tools used by the default kernel.
CONDA_HOST_REQUIREMENTS = "hdf5 * nompi_*" zlib libopenblas "fftw * nompi_*"
CONDA_RUN_REQUIREMENTS = hdf5 zlib libopenblas "fftw * nompi_*"
