# Module packaging options

Each source module or plugin may have a manually maintained `conda.mk` in its
directory. Omit the file when no packaging declarations are needed. The recipe
generator never rewrites these inputs. Run `python extensions/conda_recipe.py`
from the repository root after editing them; the build scripts do this too.
Recipes, variant YAML files and the selected-regression JSON are generated and
are not committed. Platform defaults for the base/extension recipes remain in
their `conda_build_config.yaml.in` templates.

This is a small declaration language, not an executable Makefile. Use one
assignment per line. `=` replaces a list, `+=` appends, and `#` introduces a
comment outside quotes. Quote a dependency containing spaces, and quote each
complete test command. Unknown settings and malformed values are errors.
There are no Make includes, substitutions, continuations or shell expansions.

```make
# Options: see conda/module-options.md (repository root).
CONDA_HOST_REQUIREMENTS = libboost-devel
CONDA_RUN_REQUIREMENTS = libboost
CONDA_TEST_REGTESTS = drr/rt-checkpoint-drr
```

## Requirements and package selection

| Setting | Default | Meaning |
| --- | --- | --- |
| `CONDA_BUILD_REQUIREMENTS` | Empty | Tools executed on the build machine, in addition to the backend's standard compiler/make/Python tools. |
| `CONDA_HOST_REQUIREMENTS` | Empty | Target libraries and headers used to compile this module. |
| `CONDA_RUN_REQUIREMENTS` | Empty | Runtime packages. The generator adds compatible ABI pins and preserves explicit version/build restrictions. |
| `CONDA_TEST_REQUIREMENTS` | Empty | Extra packages installed for this output's conda-build tests. |
| `CONDA_UNAVAILABLE_LIBRARIES` | Empty | Informational list of libraries without a verified conda-forge provider. Does not exclude the module. |
| `CONDA_PLATFORMS` | Unrestricted | Allowed conda subdirectories, such as `linux-64 linux-aarch64 osx-64 osx-arm64`. Intersects the backend's supported platforms. |
| `CONDA_DUMMY_PACKAGE` | `no` | `yes` preserves `plumed-ext-NAME` as an empty compatibility package when the module is built into the kernel. |

Requirements accept trailing conda selectors. A selector applies to all entries
on that line; MatchSpecs containing spaces must be quoted:

```make
CONDA_HOST_REQUIREMENTS = "hdf5 * nompi_*"
CONDA_RUN_REQUIREMENTS = hdf5
CONDA_HOST_REQUIREMENTS += libpython  # [py>=315]
CONDA_RUN_REQUIREMENTS += libpython  # [py>=315]
```

Hard module dependencies come from native `USE` declarations; do not repeat
them here. The base package inherits declarations from `always`/`default-on`
modules and their hard dependency closure. A platform restriction on a hard
dependency also restricts packages that require it. Libraries belonging only
to default-off modules stay out of the base. Licenses and notices continue to
come from module files. Shortcut/action dependencies are not inferred.

## Variants and permitted combinations

`CONDA_VARIANT_NAME` declares a conda-build variant dimension. The part after
`CONDA_VARIANT_` is the exact conda key, including its case. For example:

```make
CONDA_VARIANT_python = 3.11 3.12 3.13 3.14 3.15
CONDA_VARIANT_numpy = 2
CONDA_VARIANT_CONDA_BUILD_SYSROOT = /opt/MacOSX10.13.sdk  # [osx and x86_64]
CONDA_VARIANT_CONDA_BUILD_SYSROOT += /opt/MacOSX11.0.sdk  # [osx and arm64]
```

All values are strings, including version numbers. Variant value lists support
selectors, just like requirements. Dimensions, their forced hash keys and zipped
groups are included only where at least one value is selected.
Different dimensions normally form a Cartesian
product. `CONDA_ZIP_KEYS` declares comma-separated groups of dimensions that must
advance together. For example, these declarations allow only two SDK/CUDA pairs:

```make
CONDA_VARIANT_nvhpc_version = 25.11 25.11
CONDA_VARIANT_openacc_cuda_version = 12.9 13.0
CONDA_ZIP_KEYS = nvhpc_version,openacc_cuda_version
```

This illustrates syntax; a real SDK must contain every CUDA version requested.
Zipped lists must have equal lengths and matching selectors. Multiple groups
can be listed as separate words. A module-owned dimension overrides the same
dimension in its platform template. Conflicting dimensions in a shared recipe
are errors. Native extension outputs share one recipe/build environment;
their dimensions are combined, rather than providing independent build systems
for each output. Python/CUDA/OpenACC plugins have separate recipes. Variant keys
participate in build hashes to distinguish alternative artifacts. All variants
continue to pin the exact main kernel; no per-extension MPI matrix is needed.

## Build options and feature policy

| Setting | Default | Meaning |
| --- | --- | --- |
| `CONDA_CONFIGURE_FLAGS` | Empty | Additional arguments for the kernel's configure command. |
| `CONDA_CPPFLAGS` | Empty | Additional preprocessor flags. |
| `CONDA_CXXFLAGS` | Empty | Additional C++ compilation flags. |
| `CONDA_LDFLAGS` | Empty | Additional link flags. |
| `CONDA_MAKE_ARGS` | Empty | Extra Make arguments for the kernel or Make-based plugin backend. |
| `CONDA_CMAKE_ARGS` | Empty | CMake arguments for the pycv/scikit-build backend. |
| `CONDA_REQUIRED_FEATURES` | Empty | Features that must be available; a failed probe/check fails the build. |
| `CONDA_OPTIONAL_FEATURES` | Empty | Additional explicitly requested optional feature probes. Source guards are still probed automatically. Absence leaves the module enabled. |
| `CONDA_DISABLED_FEATURES` | Empty | Features to disable deliberately and report as absent. |

Flags are lists of arguments; quote an argument containing spaces. For build
flags and configure/Make arguments, `$PREFIX`, `${PREFIX}` and other existing
build-environment variables are expanded once by the build helper. `$$` means
a literal dollar. Missing variables are errors. Values are never evaluated as
shell commands. Test commands are different: conda-build executes them normally
in its test shell, so shell variables and conda/Jinja expressions may be used.

Native extension builds apply CPP/CXX/link flags per module and preserve the
kernel's ABI flags. Those builds use a Python compiler driver, so use the flag
and feature settings rather than configure/Make/CMake options. Kernel modules
share one configure/compiler context: their settings are combined globally.
The CUDA/OpenACC configuration scripts currently have no configurable switches;
use flags or `CONDA_MAKE_ARGS` for those plugins. CUDA flags must be accepted by
NVCC; backend-specific variables such as `NVCCCFLAGS` can be supplied through
`CONDA_MAKE_ARGS`. OpenACC link flags extend its NVIDIA shared-link command.
Avoid overriding ABI/module-selection options or disabling kernel-wide required
facilities. Changing a module between kernel and extension remains a
`module.type`/native build decision, not a packaging override.

Feature names are uppercase suffixes of `__PLUMED_HAS_*`, for example
`BOOST_SERIALIZATION`, with no `boost` alias. Native optional-library probes
currently support `BOOST_SERIALIZATION`, `FFTW`, `LIBTORCH`, `LIBMETATOMIC`,
`LIBHDF5`, `LIBITENSOR` and `KENREF`. Kernel requests use the corresponding
configure options; plugin required/disabled declarations are checked against
the installed module's exported feature vector. A required/optional feature
cannot also be disabled. Unknown native probes are errors.

The proof-of-concept default remains **all libraries optional**: source guards
are probed as before, and an unavailable optional dependency does not disable a
module. Setting `CONDA_REQUIRED_FEATURES` explicitly opts into a stricter policy;
this work does not automatically reclassify existing modules.

## Tests and hardware requirements

| Setting | Default | Meaning |
| --- | --- | --- |
| `CONDA_TEST_COMMANDS` | Empty | Quoted shell commands appended to this package's tests. |
| `CONDA_TEST_SOURCE_FILES` | Empty | Repository-relative files/directories copied into its conda-build test directory. |
| `CONDA_TEST_REGTESTS` | Empty | Existing `regtest/` fixtures selected for the installed-package integration stage. Paths are relative to `regtest/`. |
| `CONDA_TEST_PACKAGES` | Empty | Additional package names needed by those selected regressions. The owning module package is included automatically. |
| `CONDA_TEST_FORCE_REQUIREMENTS` | `no` | `yes` bypasses eligibility checks for that explicitly curated selection. Calculations/comparisons still run, and skipped tests fail package validation. |
| `CONDA_TEST_HARDWARE` | `none` | `gpu` marks custom commands and selected regressions as requiring GPU hardware. |

```make
CONDA_TEST_REQUIREMENTS = python
CONDA_TEST_COMMANDS = 'python extensions/tests/openacc.py --prefix "$PREFIX" --device host'
CONDA_TEST_SOURCE_FILES = extensions/tests/openacc.py regtest/scripts
```

Supply all fixtures required by a command; the short example above is not the
complete OpenACC source list. See the actual plugin `conda.mk` files for complete
working declarations. Loading/presence tests remain unconditional. GPU-marked
custom commands are included only when recipes are generated with
`PLUMED_CONDA_GPU_TESTS=yes`; GPU-marked selected regressions require the
installed regression runner's `--gpu` option. This is an explicit hardware-test
selection, not a claim that ordinary CPU CI validates GPU execution.

Selected numerical fixtures are collected into generated
`conda/plumed-extensions/regtests.json` and exercised by
`extensions/tests/regressions.py` during `make -C conda test`. An additional
runtime package needed by that stage belongs in `CONDA_TEST_PACKAGES`;
`CONDA_TEST_REQUIREMENTS` belongs to the separate conda-build test environment.
Each fixture's normal `plumed_needs="module:feature"` requirements are checked
before reading its input. Prefer normal checks; forcing is an explicit exception.

After editing, run the configuration tests and regenerate:

```sh
python conda/tests/module_options.py
python extensions/conda_recipe.py
```
