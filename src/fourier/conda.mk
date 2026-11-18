# Options: see conda/module-options.md (repository root).
CONDA_HOST_REQUIREMENTS = "fftw * nompi_*"
CONDA_RUN_REQUIREMENTS = "fftw * nompi_*"

CONDA_TEST_REGTESTS = fourier/rt-fftw
CONDA_TEST_PACKAGES = plumed-ext-symfunc
