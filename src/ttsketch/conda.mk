# Build against serial HDF5, but allow ABI-compatible MPI variants at runtime.
CONDA_HOST_REQUIREMENTS = "hdf5 * nompi_*"
CONDA_RUN_REQUIREMENTS = hdf5
# Keep the module enabled without ITensor while compulsory libraries are reviewed.
CONDA_UNAVAILABLE_LIBRARIES = ITensor
