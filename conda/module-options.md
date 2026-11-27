# Module declarations and packaging options

Each source module has a `module.yml`; plugins use `plugin.yml`. These files
contain descriptive metadata and an optional `conda:` mapping. Recipes, variant
YAML files and selected-regression JSON are generated and are not committed.
After editing a declaration, run:

```sh
python -m pip install -r extensions/requirements.txt
python conda/tests/module_options.py
python extensions/conda_recipe.py
```

Python metadata and extension tools require PyYAML. Ordinary configure/Make and
CMake builds do not require Python or a YAML library.

## Native build fields

Source modules declare their selection policy and direct native dependencies:

```yaml
name: drr
type: default-off
use: [core, tools, bias, cltools, config]
```

`type` must be `always`, `default-on` or `default-off`. `use` replaces the former
Makefile `USE=` declaration, with the same dependency order and meaning. These
are source dependencies, not action/shortcut requirements. Keep `config` where
needed; it is generated build support rather than a selectable module. Cycles
among infrastructure modules are permitted. Use `use: []` for no dependencies.
Makefiles still contain executable build rules and module-specific flags.

Two optional fields say how a module is compiled and where it may be placed:

```yaml
toolchain: cuda
placement: any
```

`toolchain` is `cxx` (the default: every source is compiled by the C++ compiler)
or `cuda`: the `.cu` files are compiled by nvcc, or by hipcc for AMD GPUs when
PLUMED is configured with `--enable-hip`, and the `.cpp` files by the C++
compiler. Either way the library is linked by the C++ compiler, so such a module
can go into the kernel. `placement` is `any` (the default: into the kernel with
`--enable-modules`, or as an extension against an installed base) or `extension`
(only as an extension). `--enable-modules=all` enables neither modules with a
toolchain other than `cxx` nor `placement: extension` modules; a CUDA module is
enabled only by name (`--enable-modules=+gpu_cuda_hip`), and only then does
configure look for nvcc or hipcc. The CMake build never compiles them.

A module that uses only the headers of another module lists it in `include`:

```yaml
use: [core, tools, colvar]
include: [gpu_cuda_hip]
```

`include` makes the headers of the listed modules visible to the sources of the
module, and nothing else: unlike `use`, an included module is not enabled, linked or
required at runtime. So `--enable-modules=+gpu_voronoi` puts only gpu_voronoi into
the kernel, `libplumed-ext-gpu_voronoi` does not link `libplumed-ext-gpu_cuda_hip`, and
its conda package does not depend on that of gpu_cuda_hip. This works only if the
included headers are self-contained: everything in them must be inline, templates
(including CUDA kernels) or macros, so that no object of the included module is ever
needed, and they can include only headers of modules that the including module
uses. The headers of gpu_cuda_hip used in this way (`cuda_to_hip.cuh` and
`cudaHelpers.cuh`) follow this rule.

For portable readers, the build fields have a deliberately restricted syntax:
exactly one top-level line `type: VALUE` and one `use: [name, other-name]`, at
most one `toolchain: VALUE`, one `placement: VALUE` and one `include: [name]`, with
one space after the colon and comma, no quotes, trailing comments, aliases or
multiline forms. Dependency identifiers match `[a-z][a-z0-9_-]*`. The shared
POSIX reader is `src/maketools/module-header`; Python and CMake enforce the same
grammar. The remaining metadata uses ordinary YAML, without aliases or merge
keys; duplicate mapping keys are errors. Shell/Make expressions are never evaluated.

Infrastructure modules may contain only `name`, `type` and `use`. Adding their
build declarations does not create user-manual pages; the manual still selects
module pages from action metadata. Documented modules retain `authors`,
`description`, `dois`, `tags` and their `module.md`. Plugins retain their existing
documentation fields and separate build backends; they do not use native module
selection `type`/`use` fields. License information comes from license files.

## Conda requirements

```yaml
conda:
  requirements:
    build: [cmake, ninja]
    host:
      - 'hdf5 * nompi_*'
      - spec: libpython
        selector: 'py>=315'
    run: [hdf5]
    test: [python]
  platforms: [linux-64, linux-aarch64, osx-64, osx-arm64]
  dummy_package: false
  unavailable_libraries: [example-library]
```

`requirements.build` lists build-machine tools, `host` target headers/libraries,
`run` runtime packages, and `test` extra packages for conda-build tests. Each
entry is a string or a `{spec: STRING, selector: STRING}` mapping. Selectors are
explicit data, not YAML comments. Compatible ABI pins are added to runtime
requirements while explicit constraints are retained.

The omitted `conda` section means empty declarations, unrestricted platforms
and no dummy output. `platforms` restricts allowed conda subdirectories and
intersects backend support. `dummy_package: true` preserves an empty
`plumed-ext-NAME` compatibility output for a module built into the kernel.
`unavailable_libraries` is informational and does not disable the module.

Hard package dependencies come from `use`, without duplication in `conda`.
The base package inherits requirements of `always`/`default-on` modules and
their native dependency closure. Restrictions propagate through hard dependencies.
Libraries belonging only to default-off modules stay out of the base package.

## Variants

```yaml
conda:
  variants:
    python: ['3.11', '3.12', '3.13', '3.14', '3.15']
    CONDA_BUILD_SYSROOT:
      - spec: /opt/MacOSX11.0.sdk
        selector: 'osx and arm64'
    nvhpc_version: ['25.11', '25.11']
    openacc_cuda_version: ['12.9', '13.0']
  zip_keys: ['nvhpc_version,openacc_cuda_version']
```

All variant values must be strings: quote versions and numeric values. Variant
entries accept the same spec/selector mappings as requirements. Keys preserve
case. Dimensions normally form a Cartesian product; each comma-separated
`zip_keys` group advances together and must have equal lengths and selectors.
The example illustrates syntax, not verified SDK/CUDA combinations.

Module dimensions override platform template dimensions. Conflicting dimensions
in a shared recipe are errors. Native extension outputs share one build
context; Python/CUDA/OpenACC plugins have separate recipes. Variant keys take
part in build hashes; every extension pins the exact main kernel build.
Platform defaults remain in `conda_build_config.yaml.in` templates.

## Build options and features

All the following `conda` fields default to empty lists of strings:

| Field | Meaning |
| --- | --- |
| `configure_flags` | Extra kernel configure arguments. |
| `cppflags`, `cxxflags`, `ldflags` | Extra compiler/preprocessor/linker arguments. |
| `make_args` | Kernel or Make-based plugin arguments. |
| `cmake_args` | pycv/scikit-build arguments. |
| `required_features` | Missing feature fails the build. |
| `optional_features` | Explicit additional optional probes. |
| `disabled_features` | Deliberately disable and report features absent. |

```yaml
conda:
  cppflags: ['-I${PREFIX}/include path']
  required_features: [BOOST_SERIALIZATION]
```

Flags are argument lists. Build-environment variables expand once; `$$` means a
literal dollar. Missing variables are errors. Values are never shell-evaluated.
Kernel settings combine globally. Native extensions apply flags per module,
retain kernel ABI flags and use a Python compiler driver, so configure/Make/CMake
arguments are rejected for that backend. pycv accepts `cmake_args`; CUDA/OpenACC
accept `make_args`. Kernel configure arguments cannot override installation
layout, module selection or MPI policy.

Feature names are uppercase `__PLUMED_HAS_` suffixes, e.g.
`BOOST_SERIALIZATION`, without a `boost` alias. Native probes support
`BOOST_SERIALIZATION`, `FFTW`, `LIBTORCH`, `LIBMETATOMIC`, `LIBHDF5`, `LIBITENSOR`
and `KENREF`. Unknown probes and contradictory feature policies are errors.
Libraries remain optional by default: absent optional libraries do not disable
modules. Explicit required features opt into stricter checks. Runtime inquiries
continue to report compiled feature provenance.

## Tests

```yaml
conda:
  requirements:
    test: [python]
  tests:
    regtests: [drr/rt-checkpoint-drr]
    packages: []
    commands: ['python extensions/tests/example.py --prefix "$PREFIX"']
    source_files: [extensions/tests/example.py]
    force_requirements: false
    hardware: none
```

`commands` contains shell commands executed by conda-build. `source_files`
contains repository-relative fixtures copied to its test directory. `regtests`
selects installed integration fixtures relative to `regtest/`; `packages`
adds packages needed by that integration stage, with the owning module package
included automatically. `requirements.test` applies to the separate conda-build
test environment. Lists default empty.

`force_requirements` defaults false: true bypasses eligibility checks for a
curated selection, while calculations/comparisons still run and skips fail.
`hardware` defaults `none`; `gpu` requires explicit GPU selection. Recipes include
GPU custom commands only with `PLUMED_CONDA_GPU_TESTS=yes`; the installed
regression runner requires `--gpu` for GPU selections. Presence/loading tests
remain unconditional. Ordinary CPU CI does not establish GPU execution.

Selected fixtures go into generated `conda/plumed-extensions/regtests.json` and
run through `extensions/tests/regressions.py` during `make -C conda test`.
Keep ordinary module/feature eligibility checks unless forcing is intentional.
