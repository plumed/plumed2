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
#ifndef __PLUMED_cuda_helpers_cuh
#define __PLUMED_cuda_helpers_cuh

#include "plumed/tools/Tensor.h"
#include "plumed/tools/Vector.h"
#include "plumed/tools/View.h"
#include "plumed/tools/Exception.h"
#include <cub/block/block_reduce.cuh>
#include <thrust/device_vector.h>
#include <algorithm>
#include <string>
#include <vector>
namespace CUDAHELPERS {

/// @brief throws a plumed exception if a CUDA runtime call failed
inline void cudaCheck (cudaError_t const err, const char *what) {
  if (err != cudaSuccess) {
    plumed_merror (std::string (what) + ": " + cudaGetErrorString (err));
  }
}

/// @brief a minimal grow-only, uninitialized buffer in device memory
///
/// Unlike thrust::device_vector this never launches initialization kernels and
/// never works on the legacy default stream: the only (implicitly
/// device-synchronizing) calls are cudaMalloc/cudaFree, that happen only when
/// the buffer has to grow.
/// The content is NOT preserved when the buffer grows.
template <typename T>
class deviceBuffer {
  T *ptr_ = nullptr;
  size_t size_ = 0;
  size_t capacity_ = 0;
public:
  deviceBuffer() = default;
  deviceBuffer (const deviceBuffer &) = delete;
  deviceBuffer &operator= (const deviceBuffer &) = delete;
  ~deviceBuffer() {
    if (ptr_) {
      cudaFree (ptr_);
    }
  }
  /// @param slack the capacity is set to n*slack when the buffer has to grow
  void resize (size_t const n, double const slack = 1.0) {
    if (n > capacity_) {
      if (ptr_) {
        cudaCheck (cudaFree (ptr_), "cudaFree");
        ptr_ = nullptr;
        capacity_ = 0;
      }
      const size_t newCapacity = std::max (n, static_cast<size_t> (n * slack));
      cudaCheck (cudaMalloc (&ptr_, newCapacity * sizeof (T)), "cudaMalloc");
      capacity_ = newCapacity;
    }
    size_ = n;
  }
  T *data() {
    return ptr_;
  }
  const T *data() const {
    return ptr_;
  }
  size_t size() const {
    return size_;
  }
};

/// @brief a minimal grow-only, uninitialized buffer in pinned (page-locked) host memory
///
/// Pinned memory makes cudaMemcpyAsync really asynchronous and faster.
/// The content is NOT preserved when the buffer grows.
template <typename T>
class pinnedBuffer {
  T *ptr_ = nullptr;
  size_t size_ = 0;
  size_t capacity_ = 0;
public:
  pinnedBuffer() = default;
  pinnedBuffer (const pinnedBuffer &) = delete;
  pinnedBuffer &operator= (const pinnedBuffer &) = delete;
  ~pinnedBuffer() {
    if (ptr_) {
      cudaFreeHost (ptr_);
    }
  }
  /// @param slack the capacity is set to n*slack when the buffer has to grow
  void resize (size_t const n, double const slack = 1.0) {
    if (n > capacity_) {
      if (ptr_) {
        cudaCheck (cudaFreeHost (ptr_), "cudaFreeHost");
        ptr_ = nullptr;
        capacity_ = 0;
      }
      const size_t newCapacity = std::max (n, static_cast<size_t> (n * slack));
      cudaCheck (cudaMallocHost (&ptr_, newCapacity * sizeof (T)), "cudaMallocHost");
      capacity_ = newCapacity;
    }
    size_ = n;
  }
  T *data() {
    return ptr_;
  }
  const T *data() const {
    return ptr_;
  }
  size_t size() const {
    return size_;
  }
  T *begin() {
    return ptr_;
  }
  T *end() {
    return ptr_ + size_;
  }
  T &operator[] (size_t const i) {
    return ptr_[i];
  }
  const T &operator[] (size_t const i) const {
    return ptr_[i];
  }
};

/// The dynamic shared memory of the kernel, seen as an array of T
///
/// For kernels that use the dynamic shared memory for a single purpose;
/// kernels that split it in more arrays should use sharedArena
template <typename T>
__device__ __forceinline__ T *dynamicSharedMemory() {
  extern __shared__ __align__ (16) unsigned char memory[];
  return reinterpret_cast<T *> (memory);
}

///A simple struct to use the allocated shared memory in multiple ways
///
///Usage: declare this at the begininnig of your kernel (a single instance per kernel)
///and then "allocate" memory by calling the `get_shared_memory method`
class sharedArena {
  //TODO: this costs a register, it may be a better idea to use an external register and then override it for another usage
  unsigned allocated=0;
public:
  template <typename T>
  __device__ T *get_shared_memory(unsigned const size) {
    // the caller is responsible for the alignment of the subsequent arrays
    auto ptr = reinterpret_cast<T *> (dynamicSharedMemory<unsigned char>() + allocated);
    allocated+=size * sizeof(T);
    return ptr;
  }
};

/// @brief reduces N values per thread over the block and writes the N block
/// totals in partials[k*nBlocks + blockId]
///
/// Deterministic (no atomics). All the threads of the block must call this,
/// blockDim.x must be a multiple of 32 and smem must hold at least
/// N*(blockDim.x/32) elements. smem can alias other shared memory as long as
/// the caller has synchronized the block after its last use.
template <unsigned N, typename T>
__device__ __forceinline__ void blockReduceToPartials (const T (&vals)[N],
    T *smem,
    T *partials,
    const unsigned blockId,
    const unsigned nBlocks) {
  const unsigned lane = threadIdx.x % 32u;
  const unsigned warp = threadIdx.x / 32u;
  const unsigned nwarps = blockDim.x / 32u;
#pragma unroll
  for (unsigned k = 0; k < N; ++k) {
    T v = vals[k];
#pragma unroll
    for (unsigned offset = 16; offset > 0; offset /= 2) {
      v += __shfl_down_sync (0xffffffffu, v, offset);
    }
    if (lane == 0) {
      smem[warp * N + k] = v;
    }
  }
  __syncthreads();
  if (threadIdx.x < N) {
    T acc = T (0);
    for (unsigned w = 0; w < nwarps; ++w) {
      acc += smem[w * N + threadIdx.x];
    }
    partials[threadIdx.x * nBlocks + blockId] = acc;
  }
}

/// @brief shared memory needed by blockReduceToPartials
template <unsigned N, typename T>
constexpr size_t blockReduceMemory (unsigned const threads) {
  return N * ((threads + 31) / 32) * sizeof (T);
}

/// @brief final step of blockReduceToPartials: to be launched with one block
/// per accumulated quantity, out[blockIdx.x] = sum_i partials[blockIdx.x*nPartials+i]
constexpr unsigned reducePartialsThreads = 256;
template <typename T, unsigned THREADS = reducePartialsThreads>
__global__ void reducePartials (const T *partials,
                                const unsigned nPartials,
                                T *out) {
  using BlockReduceT = cub::BlockReduce<T, THREADS>;
  __shared__ typename BlockReduceT::TempStorage temp_storage;
  const T *in = partials + blockIdx.x * nPartials;
  T acc = T (0);
  for (unsigned i = threadIdx.x; i < nPartials; i += THREADS) {
    acc += in[i];
  }
  const T total = BlockReduceT (temp_storage).Sum (acc);
  if (threadIdx.x == 0) {
    out[blockIdx.x] = total;
  }
}

} // namespace CUDAHELPERS
#endif //__PLUMED_cuda_helpers_cuh
