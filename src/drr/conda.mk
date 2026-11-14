# Options: see conda/module-options.md (repository root).
# Boost serialization is available on conda-forge; it is optional for this module.
CONDA_HOST_REQUIREMENTS = libboost-devel
CONDA_RUN_REQUIREMENTS = libboost

CONDA_TEST_REGTESTS = drr/rt-checkpoint-drr
