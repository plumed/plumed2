# Options: see conda/module-options.md (repository root).
# Driver trajectory formats; chemfiles has no run_exports and changes SONAME
# with its minor version. 0.10.3 also supports linux-aarch64.
CONDA_HOST_REQUIREMENTS = "chemfiles-lib 0.10.*"
CONDA_RUN_REQUIREMENTS = "chemfiles-lib >=0.10.3,<0.11"
CONDA_CONFIGURE_FLAGS = --enable-libchemfiles
