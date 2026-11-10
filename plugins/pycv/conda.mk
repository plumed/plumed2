# Python variants are selected by conda/plumed-pycv/conda_build_config.yaml.
CONDA_HOST_REQUIREMENTS = python "pybind11 >=3.0" numpy
CONDA_RUN_REQUIREMENTS = python numpy
# Python 3.15 splits its embedding library into a separate package.
CONDA_HOST_REQUIREMENTS += libpython  # [py>=315]
CONDA_RUN_REQUIREMENTS += libpython  # [py>=315]
