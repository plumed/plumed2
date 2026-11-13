# Options: see conda/module-options.md (repository root).
# Python variants are selected by conda/plumed-pycv/conda_build_config.yaml.
CONDA_HOST_REQUIREMENTS = python "pybind11 >=3.0" numpy
CONDA_RUN_REQUIREMENTS = python numpy
# Python 3.15 splits its embedding library into a separate package.
CONDA_HOST_REQUIREMENTS += libpython  # [py>=315]
CONDA_RUN_REQUIREMENTS += libpython  # [py>=315]

CONDA_PLATFORMS = linux-64 linux-aarch64 osx-64 osx-arm64
CONDA_BUILD_REQUIREMENTS = cmake ninja
CONDA_VARIANT_python = 3.11 3.12 3.13 3.14 3.15
CONDA_VARIANT_numpy = 2
CONDA_VARIANT_CONDA_BUILD_SYSROOT = /opt/MacOSX10.13.sdk  # [osx and x86_64]
CONDA_VARIANT_CONDA_BUILD_SYSROOT += /opt/MacOSX11.0.sdk  # [osx and arm64]
CONDA_VARIANT_c_stdlib = sysroot  # [linux]
CONDA_VARIANT_c_stdlib_version = 2.28  # [linux]
CONDA_TEST_REQUIREMENTS = pytest py-plumed
CONDA_TEST_SOURCE_FILES = extensions/tests/pycv.py plugins/pycv/pythontests plugins/pycv/regtest regtest/scripts
CONDA_TEST_COMMANDS = 'python extensions/tests/pycv.py --prefix "$PREFIX"'
