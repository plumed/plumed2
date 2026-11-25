/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2024 Daniele Rapetti, The plumed team

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   cudaOnPlumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   cudaOnPlumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_cuda_to_hip_cuh
#define __PLUMED_cuda_to_hip_cuh

// The CUDA runtime names used by the plugin (and by ../openaccPTM/gpucheck),
// mapped to HIP when the plugin is compiled for AMD GPUs (USE_HIP=1 in the
// Makefile); with nvcc this only includes the CUDA runtime.
// HEADER-ONLY, as cudaHelpers.cuh: other modules include it without linking gpu_cuda_hip.

#ifdef __HIP__

// the __shfl_*_sync functions are disabled by default from ROCm 6.2 to 6.4
// (always enabled from ROCm 7.0)
#ifndef HIP_ENABLE_WARP_SYNC_BUILTINS
#define HIP_ENABLE_WARP_SYNC_BUILTINS
#endif
#include <hip/hip_runtime.h>
#include <hipcub/hipcub.hpp>

#define cudaDeviceProp                   hipDeviceProp_t
#define cudaError_t                      hipError_t
#define cudaEvent_t                      hipEvent_t
#define cudaFuncAttributes               hipFuncAttributes
#define cudaStream_t                     hipStream_t

#define cudaEventDisableTiming           hipEventDisableTiming
#define cudaMemcpyDeviceToHost           hipMemcpyDeviceToHost
#define cudaMemcpyHostToDevice           hipMemcpyHostToDevice
#define cudaStreamNonBlocking            hipStreamNonBlocking
#define cudaSuccess                      hipSuccess

#define cudaDeviceGetStreamPriorityRange hipDeviceGetStreamPriorityRange
#define cudaDeviceSynchronize            hipDeviceSynchronize
#define cudaEventCreateWithFlags         hipEventCreateWithFlags
#define cudaEventDestroy                 hipEventDestroy
#define cudaEventRecord                  hipEventRecord
#define cudaFree                         hipFree
#define cudaFreeHost                     hipHostFree
#define cudaGetDevice                    hipGetDevice
#define cudaGetDeviceCount               hipGetDeviceCount
#define cudaGetDeviceProperties          hipGetDeviceProperties
#define cudaGetErrorString               hipGetErrorString
#define cudaGetLastError                 hipGetLastError
#define cudaMalloc                       hipMalloc
#define cudaMallocHost                   hipHostMalloc
#define cudaMemcpy                       hipMemcpy
#define cudaMemcpyAsync                  hipMemcpyAsync
#define cudaMemset                       hipMemset
#define cudaMemsetAsync                  hipMemsetAsync
#define cudaSetDevice                    hipSetDevice
#define cudaStreamCreateWithFlags        hipStreamCreateWithFlags
#define cudaStreamCreateWithPriority     hipStreamCreateWithPriority
#define cudaStreamDestroy                hipStreamDestroy
#define cudaStreamSynchronize            hipStreamSynchronize
#define cudaStreamWaitEvent              hipStreamWaitEvent

// hipFuncGetAttributes takes the kernel as a const void*, CUDA also accepts
// a pointer to the kernel
template <typename Kernel>
inline hipError_t cudaFuncGetAttributes (hipFuncAttributes *attr, Kernel *kernel) {
  return hipFuncGetAttributes (attr, reinterpret_cast<const void *> (kernel));
}

namespace cub = hipcub;

#else

#include <cuda_runtime.h>

#endif

namespace CUDAHELPERS {
/// @brief the mask of the 32 lanes of the warp of the calling thread, for the
/// warp intrinsics: the kernels work with warps of 32 threads, and a 64-wide
/// AMD wavefront holds two of them
__device__ __forceinline__ auto warpMask() {
#ifdef __HIP__
  return 0xffffffffull << (__lane_id() & ~31u);
#else
  return 0xffffffffu;
#endif
}
} // namespace CUDAHELPERS

#endif //__PLUMED_cuda_to_hip_cuh
