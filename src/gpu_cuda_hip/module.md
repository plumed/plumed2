This module contains the actions that run on a GPU with CUDA (NVIDIA GPUs) or HIP (AMD GPUs):
`CUDACOORDINATION` and `CUDACOORDINATIONFLOAT`, a version of
[COORDINATION](COORDINATION.md) that runs on the GPU, in double and single precision.
The same CUDA sources are compiled with nvcc for NVIDIA GPUs and with hipcc for AMD GPUs.

These actions were developed in the lesson of the [PLUMED school](https://plumed-school.github.io/lessons/23/004/data/NAVIGATION.html)
that presents them step by step. Up to PLUMED 2.10 they were in the `cudaCoord` plugin and had to be loaded with
`LOAD FILE=CudaCoordination.so`: this is no longer needed, and such a line must be removed from the input.

## Installation

This module is not installed by default, and `--enable-modules=all` does not enable it: it has to be named.
On a machine with the CUDA toolkit (nvcc), add `--enable-modules=+gpu_cuda_hip` to the `./configure` command:

```bash
./configure --enable-modules=+gpu_cuda_hip
```

configure looks for nvcc (set `NVCC` to choose another one) only when this module is enabled. The `.cu` files are
compiled by nvcc, with the flags of the C++ compiler, and the kernel is linked by the C++ compiler as usual.
Useful options:

- `NVCCFLAGS`: additional flags for nvcc, e.g. `NVCCFLAGS=-arch=sm_80` to compile for a given GPU (by default nvcc
  compiles for its default architecture, and the code is compiled at runtime for the GPU in use);
- `NVCC_CCBIN`: the host compiler of nvcc, by default the C++ compiler of PLUMED (ccache is removed);
- `--with-cudart=static|shared`: the CUDA runtime linked into the kernel. With the static runtime (the default) the kernel
  can be loaded also on machines without the CUDA toolkit or a GPU, where these actions stop with an error;
- `CUDA_LIBDIR`: the directory of the CUDA runtime library, if configure does not find it.

For AMD GPUs, with ROCm 6.2 or later (hipcc, rocThrust and hipCUB) installed, add `--enable-hip`:

```bash
./configure --enable-modules=+gpu_cuda_hip --enable-hip HIP_ARCHITECTURES="gfx90a gfx942"
```

`HIP_ARCHITECTURES` lists the GPU targets (by default the GPUs of the machine that compiles), `HIPCC` chooses the compiler.
The HIP runtime is a shared library: the kernel then needs ROCm to be loaded.

With conda, this module is installed separately, in a package for each CUDA version, on top of the same `plumed` package.

## Tests

The tests are in `regtest/gpu_cuda_hip`. The tests in `nogpu` (input errors, the error when there is no GPU) run everywhere
the module is compiled; the others need a GPU, and are reported as not applicable where `plumed cuda-devices` finds none.
On a node with GPUs, `plugins/openaccPTM/gpucheck` checks that the actions leave the GPU of the MD code alone.

## Code shared with other modules

`cuda_to_hip.cuh` (the CUDA names mapped to HIP) and `cudaHelpers.cuh` (GPU choice, device guard, buffers, deterministic
reductions) are also used by other CUDA/HIP modules (e.g. [gpu_voronoi](module_gpu_voronoi.md)), which list
`include: [gpu_cuda_hip]` in their `module.yml`: they use these headers without enabling or linking this module, so they
must stay header-only (inline functions, templates and kernels, macros).

## Limitations

`CUDACOORDINATION` and `CUDACOORDINATIONFLOAT` work as [COORDINATION](COORDINATION.md),
except that:

- the minimum image convention is applied in scaled coordinates: it is exact for an orthorhombic box in any orientation
  (e.g. the box rotated by [FIT_TO_TEMPLATE](FIT_TO_TEMPLATE.md)); for a skewed (triclinic) box it is exact only up to half
  the smallest height of the box, so `D_MAX` must be smaller than that (an error is raised otherwise);
- the SWITCH keyword is not supported;
- the THREADS keyword sets the maximum number of threads per block used by the kernels (rounded down to a multiple of 32,
  the cell kernels use at most 256); there is no limit on the number of atoms per cell.

## Which GPU is used

The GPU is chosen by PLUMED (`PLMD::GpuDevice`, in `src/core/GpuDevice.h`), with the same rules as the other actions that
run on a GPU: the GPU assigned to PLUMED by the MD code (with the `setGpuDeviceId` command, as the GROMACS and NAMD patches
do when GROMACS or NAMD use CUDA), otherwise the GPU the MD code is using (the one whose CUDA context is current when the
action is created), otherwise, with several GPUs, the rank of the process in its node modulo the number of GPUs, otherwise
the current one of the CUDA runtime (by default the first one; with `CUDA_VISIBLE_DEVICES` each rank can be given its own GPU).
It is written in the log with the rule that was applied (e.g. `c runs on GPU 1 (assigned by the MD code)` for an action
labelled `c`). With MPI only the first rank uses the GPU. Every step runs on that GPU, and the MD code finds its own GPU as it
left it, also when it is using another one.
With HIP (AMD GPUs) the module shares the runtime with the MD code, and the GPU the MD code is using is the current one of
the HIP runtime: the rank is not used, and without an MD code (`plumed driver`) every process uses the first GPU, unless
`HIP_VISIBLE_DEVICES` gives each rank its own GPU.

## Performance notes

These actions are designed to coexist with an MD code that keeps the GPU busy (e.g. GPU-resident GROMACS):

- all the work is enqueued on non-blocking, highest-priority streams, so it never implicitly synchronizes with the legacy
  default stream (and so with the other code's streams), and there are no device-wide synchronizations: each step waits
  only for its own work;
- device and pinned host buffers are persistent and grow-only, so there are no allocations (that synchronize the whole
  device) in the steady state;
- each step does one host-to-device copy (positions), the kernels, a deterministic two-stage reduction of the coordination
  and of the virial, and one device-to-host copy (derivatives, coordination and virial);
- the kernels use a warp per atom: the lanes split the candidate neighbours of the atom, and the atoms of a crowded cell
  are spread over several blocks, so that the GPU is fully used also with a few hundred atoms in a handful of cells;
- only the atoms that appear in the cell lists are moved to and from the GPU: with a small `GROUPA` in a large `GROUPB`
  (e.g. an ion in water) these are just the atoms within `NL_CUTOFF` of the `GROUPA` atoms. PLUMED itself still gathers the
  positions and scatters the forces of all the atoms in the groups at each step, so keeping the groups small is still the
  most effective optimization;
- the cell lists are built on the CPU, in a compact format with only the non-empty cells, and are uploaded only when the
  neighbor list is updated;
- by default (`NL_CUTOFF` equal to `D_MAX`) the lists are rebuilt at every step. With `NL_CUTOFF` larger than `D_MAX` they
  are rebuilt only when an atom moved more than half of the skin (`NL_CUTOFF-D_MAX`) since the last update (or every
  `NL_STRIDE` steps, if given), with the same result;
- always set `D_MAX`: the default cutoff of the rational switching function is very large (about 6.8 `R_0` for `NN=6 MM=12`).

Since PLUMED receives the positions and returns the forces through host memory, each step still pays for the data round
trip: for small groups [COORDINATION](COORDINATION.md) may be faster.
