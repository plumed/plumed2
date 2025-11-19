# Coordination in Cuda

This is the optimized version of the lesson that I presented in the [plumed-school](https://plumed-school.github.io/lessons/23/004/data/NAVIGATION.html) with a step-by-step approach.

`CUDACOORDINATION` and `CUDACOORDINATIONFLOAT` depend on [CCCL](https://github.com/NVIDIA/cccl) which is automatically fetched by the cuda compiler (if you use nvcc, you have access to the CCCL headers).

The files `cudaHelpers.cuh` and `cudaHelpers.cu` contains a few support functions for helping in interfacing `PLMD::Vector` and `PLMD::Tensor` with Cuda's thrust,
along with the reduction functions baked with Cuda's cub building blocks and their drivers.

### Compile



##### as developer
  
With a ready to run plumed:
  - `nvcc-MakeFile.sh`
  - `make`

In this way you do not need to recompile everithing if you change part of the sources

you may need to specify the SM of your GPU by modifying the `Makefile`, for example with:

```Makefile
NVCCCFLAGS = -dc -dlto --gpu-architecture=sm_75 
NVCCLDFLAGS = -shared -dlto --gpu-architecture=sm_75 
```
If you are using a T1000

#### as a user

Running `nvcc-mklib.sh` should be enought to get the `CudaCoordination.so`.

you may need to uncomment and modifly the lines
```bash
#compile="$compile --gpu-architecture=sm_75 "
#link_command="$link_command -shared -dlto --gpu-architecture=sm_75"
```
near the bottom of the script

### How to tests

 - `cd regtest; ln -s ../../../regtest/scripts .`
 - `make check`

## Limitations

`CUDACOORDINATION` and `CUDACOORDINATIONFLOAT` work more or less as the standard `COORDINATION`, except from:

 - work only with orthogonal pbcs or no pbcs at all
 - do not support the SWITCH keyword
 - new `THREADS` keyword to control the maximum number of threads that can be used by the kernels

## Current TODO
 
 - The GPU device needs to be explicitly selected
 - Integrate the CI
