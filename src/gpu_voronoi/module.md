This module contains the GPU backend of the reactive Voronoi actions [VORONOI_COORDINATION](VORONOI_COORDINATION.md),
[VORONOI_DISTANCE](VORONOI_DISTANCE.md) and [VORONOI_POSITION](VORONOI_POSITION.md): the soft assignment of the ASSIGNED
atoms to the CENTERS and its analytic derivatives, in double precision, on NVIDIA GPUs (CUDA) or AMD GPUs (HIP). The same
CUDA sources are compiled with nvcc for NVIDIA GPUs and with hipcc for AMD GPUs. The actions are in the colvar module,
and select the backend with their BACKEND keyword:

```plumed
q: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=1 POWER=2 BACKEND=CUDA
PRINT ARG=q FILE=COLVAR
```

`BACKEND=CUDA` (also on AMD GPUs) requires this backend, and stops with an error when it is not available, there is no GPU,
or the calculation is not supported; `BACKEND=AUTO` uses it when it can, and otherwise uses the CPU and writes the reason
in the log; `BACKEND=CPU` (the default) never uses it. The calculations supported by the GPU are those with one MPI rank,
an orthorhombic box or no PBC, with all the assignments or with NLIST; AUTO uses the CPU for triclinic boxes and with
more than one MPI rank.

## Installation

This module is not installed by default, and `--enable-modules=all` does not enable it: it has to be named. On a machine
with the CUDA toolkit (nvcc), add `--enable-modules=+gpu_voronoi` to the `./configure` command, for AMD GPUs also
`--enable-hip`:

```bash
./configure --enable-modules=+gpu_voronoi
```

The options of configure for CUDA and HIP (`NVCCFLAGS`, `--with-cudart`, `HIP_ARCHITECTURES`, ...) are described on the
page of the [gpu_cuda_hip](module_gpu_cuda_hip.md) module. This module uses the CUDA/HIP code of gpu_cuda_hip (its
headers only): it does not need that module, and `--enable-modules=+gpu_voronoi` adds only this backend to PLUMED.
With conda, this module is installed separately, in a package for each CUDA version, on top of the same `plumed` package.

## Tests

The tests are in `regtest/gpu_voronoi`. The tests in `rt-nogpu-*` (the fallback of AUTO and the error of BACKEND=CUDA
when there is no GPU) run everywhere the module is compiled; the others compare the GPU with the CPU backend, and need a
GPU and the gpu_cuda_hip module (they use `plumed cuda-devices` to find the GPU): elsewhere they are reported as not
applicable. The tests of the reactive Voronoi actions on the CPU are in `regtest/basic/rt-reactive-voronoi*`.

## Which GPU is used

The GPU is chosen by PLUMED (`PLMD::GpuDevice`, in `src/core/GpuDevice.h`), with the same rules as the other actions that
run on a GPU (see the page of [gpu_cuda_hip](module_gpu_cuda_hip.md)), when the backend is first used, and its number is
written in the log of the action. Every step runs on that GPU, and the MD code finds its own GPU as it left it.

## How it works

Compatible Voronoi actions in one PLUMED instance share one GPU context for each KAPPA. Each context keeps device buffers
proportional to the number of CENTERS and ASSIGNED atoms; NLIST adds the candidate list, uploaded only when it is rebuilt
on the CPU. The kernels go over the CENTERS for each ASSIGNED atom, with a shifted softmax for the normalization. All the
sums are done in a fixed order (no atomic additions of floating-point numbers), so the results are the same at every run.
The assignment is started on the GPU in advance, so that it can overlap with the work of other actions (see PLUMED\_LAUNCH
in [the page on actions](actions.md)), and its derivatives are computed when the action is reached.

OPES, the scalar reductions, the direct distance and position terms and their virial contributions are computed on the
CPU: VORONOI_DISTANCE therefore includes a reduction over the pairs of centers on the CPU. PLUMED receives the positions
and returns the forces through host memory, so each step pays for a round trip of the data: for small selections the CPU
backend can be faster. Choose the backend with matched numerical checks and timings of the actual workload.
