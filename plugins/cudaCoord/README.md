This is the optimized version of the lesson that I presented in the [plumed-school](https://plumed-school.github.io/lessons/23/004/data/NAVIGATION.html) with a step-by-step approach.

`CUDACOORDINATION` and `CUDACOORDINATIONFLOAT` depend on [CCCL](https://github.com/NVIDIA/cccl) which is automatically fetched by the cuda compiler (if you use nvcc, you have access to the CCCL headers).

The files `cudaHelpers.cuh` and `cudaHelpers.cu` contains a few support functions for helping in interfacing `PLMD::Vector` and `PLMD::Tensor` with Cuda's thrust,
along with the reduction functions baked with Cuda's cub building blocks and their drivers.

>[!WARNING]
>Plumed may refuse to to the calculations if the GPU does not allow for the calculations to be done.
>We are currently working on a solution for this
>
>A workaround is reducing the size of the cutoff of the neigbor list or to increase the number of threads with THREADS
>Plumed will present the user wiht an error message with the suggestion.

### Compile
  
With a ready to run plumed in your path:
  - `./configure`
  - `make`

you may need to specify the SM of your GPU by modifying the `Makefile`, for example with:

```Makefile
NVCCCFLAGS = -dc -dlto --gpu-architecture=sm_75 
NVCCLDFLAGS = -shared -dlto --gpu-architecture=sm_75 
```
If you are using a T1000


### How to tests

Assuming you are using this plugin from the plugin directory of a plumed source directory:
 - `cd regtest; ln -s ../../../regtest/scripts .`
 - `make check`

## Limitations

`CUDACOORDINATION` and `CUDACOORDINATIONFLOAT` work more or less as the standard `COORDINATION`, except from:

 - work only with orthogonal pbcs or no pbcs at all
 - do not support the SWITCH keyword
 - new `THREADS` keyword to control the maximum number of threads that can be used by the kernels

## Performance notes

The plugin is designed to coexist with an MD code that keeps the GPU busy (e.g. GPU-resident GROMACS):

 - all the work is enqueued on non-blocking, highest-priority CUDA streams, so it never implicitly synchronizes with the legacy default stream (and so with the other code's streams), and there are no device-wide synchronizations (`cudaDeviceSynchronize`): each step waits only for its own work with a single `cudaStreamSynchronize`
 - device and pinned host buffers are persistent and grow-only, so there are no `cudaMalloc`/`cudaFree` (that synchronize the whole device) in the steady state
 - each step does one host-to-device copy (positions), the kernels, a deterministic two-stage reduction of the coordination and of the virial, and one device-to-host copy (derivatives, coordination and virial)
 - the cell lists are built on the CPU and uploaded only when the neighbor list is updated: by default this happens at every step, set `NL_CUTOFF` larger than `D_MAX` and `NL_STRIDE` > 1 to reduce this cost

Since PLUMED receives the positions and returns the forces through host memory, each step still pays for the data round trip: for small groups the CPU `COORDINATION` may still be faster.

## Current TODO
 
 - The GPU device needs to be explicitly selected
 - Integrate the CI
