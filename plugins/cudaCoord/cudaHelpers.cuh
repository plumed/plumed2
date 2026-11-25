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
#include <thrust/device_vector.h>
#include <vector>
namespace CUDAHELPERS {

/// @brief a interface to help in the data I/O to the GPU
///
/// It contains specialized constructors for using the correct size with
/// PLMD::VectorGeneric or PLMD::TensorGeneric along with a std::vector
/// constructor to generate the corret data structure given the size of it.
///
/// The functions and types that use DataInterface assume that ALL the data
/// contained within the passed object will be used
struct DataInterface {
  // NOT owning pointer
  double *ptr = nullptr;
  size_t size = 0;
  DataInterface() = delete;

  // VectorGeneric is a "memory map" on an n linear array
  // &vg[0] gets the pointer to the first double in memory
  // C++ vectors align memory so we can safely use the vector variant of
  // DataInterface
  template <unsigned n>
  explicit DataInterface (PLMD::VectorGeneric<n> &vg)
    : ptr (&vg[0]), size (n) {}
  // TensorGeneric is a "memory map" on an n*m linear array
  // &tg[0][0] gets  the pointer to the first double in memory
  // C++ vectors align memory so we can safely use the vector variant of
  // DataInterface
  template <unsigned n, unsigned m>
  explicit DataInterface (PLMD::TensorGeneric<n, m> &tg)
    : ptr (&tg[0][0]), size (n * m) {}
  template <typename T>
  explicit DataInterface (std::vector<T> &vt) : DataInterface (vt[0]) {
    size *= vt.size();
  }
};

/// @brief the specialized function for getting double precision data from the
/// gpu to a PLMD container
inline void plmdDataFromGPU (thrust::device_vector<double> &dvmem,
                             DataInterface data) {
  cudaMemcpy (data.ptr,
              thrust::raw_pointer_cast (dvmem.data()),
              data.size * sizeof (double),
              cudaMemcpyDeviceToHost);
}

/// @brief the specialized asyncronous function for getting double precision
/// data from the gpu to a PLMD container
inline void plmdDataFromGPU (thrust::device_vector<double> &dvmem,
                             DataInterface data,
                             cudaStream_t stream) {
  cudaMemcpyAsync (data.ptr,
                   thrust::raw_pointer_cast (dvmem.data()),
                   data.size * sizeof (double),
                   cudaMemcpyDeviceToHost,
                   stream);
}

/// @brief the specialized function for putting double precision
/// data from a PLMD container to the gpu
inline void plmdDataToGPU (thrust::device_vector<double> &dvmem,
                           DataInterface data) {
  dvmem.resize (data.size);
  cudaMemcpy (thrust::raw_pointer_cast (dvmem.data()),
              data.ptr,
              data.size * sizeof (double),
              cudaMemcpyHostToDevice);
}

/// @brief the specialized asyncronous function for putting double precision
/// data from a PLMD container to the gpu
inline void plmdDataToGPU (thrust::device_vector<double> &dvmem,
                           DataInterface data,
                           cudaStream_t stream) {
  dvmem.resize (data.size);
  cudaMemcpyAsync (thrust::raw_pointer_cast (dvmem.data()),
                   data.ptr,
                   data.size * sizeof (double),
                   cudaMemcpyHostToDevice,
                   stream);
}

/// @brief the specialized function for getting single precision data from the
/// gpu to a PLMD container
/// @param dvmem the cuda interface to the data on the device
/// @param data the interface to the plumed-type of data
/// @param -ignored-
///
/// FAQ:why there is only a non async version?
/// the `for` loop used to convert float to double may start before the end of
/// async data transfer, so is safer to use the non async version, the stream
/// paraemter is ignored so that there is non need to write extra logic in the
/// user function tha may call the templated double version
inline void plmdDataFromGPU (thrust::device_vector<float> &dvmem,
                             DataInterface data,
                             cudaStream_t = 0) {
  std::vector<float> tempMemory (data.size);
  cudaMemcpy (tempMemory.data(),
              thrust::raw_pointer_cast (dvmem.data()),
              data.size * sizeof (float),
              cudaMemcpyDeviceToHost);
  for (auto i = 0u; i < data.size; ++i) {
    data.ptr[i] = tempMemory[i];
  }
}

// since the for loop is BEFORE the memcpy I can safely distinguish the async
// case from the non async one

/// @brief the specialized function for putting double precision
/// data from a PLMD container to the gpu
inline void plmdDataToGPU (thrust::device_vector<float> &dvmem,
                           DataInterface data) {
  dvmem.resize (data.size);
  std::vector<float> tempMemory (data.size);
  for (auto i = 0u; i < data.size; ++i) {
    tempMemory[i] = data.ptr[i];
  }
  cudaMemcpy (thrust::raw_pointer_cast (dvmem.data()),
              tempMemory.data(),
              data.size * sizeof (float),
              cudaMemcpyHostToDevice);
}

/// @brief the specialized asyncronous function for putting double precision
/// data from a PLMD container to the gpu
inline void plmdDataToGPU (thrust::device_vector<float> &dvmem,
                           DataInterface data,
                           cudaStream_t stream) {
  dvmem.resize (data.size);
  std::vector<float> tempMemory (data.size);
  for (auto i = 0u; i < data.size; ++i) {
    tempMemory[i] = data.ptr[i];
  }
  cudaMemcpyAsync (thrust::raw_pointer_cast (dvmem.data()),
                   tempMemory.data(),
                   data.size * sizeof (float),
                   cudaMemcpyHostToDevice,
                   stream);
}

// the explicit constructors of DataInterface create the need for a wrapper
/// @brief copies data to the GPU, using thrust::device_vector as interface
///
/// if you are transferring data from/to a 'standard' vector, consider using
/// thrust::host_vector<T> and the thrust syntax (`deviceV=hostV;` for moving to
/// the device or `hostV=`deviceV;` for moving from the device)
template <typename T, typename Y>
inline void plmdDataToGPU (thrust::device_vector<T> &dvmem, Y &data) {
  plmdDataToGPU (dvmem, DataInterface (data));
}

/// @brief async version of plmdDataToGPU
///
/// if you are transferring data from/to a 'standard' vector, consider using
/// thrust::host_vector<T> and the thrust syntax (`deviceV=hostV;` for moving to
/// the device or `hostV=`deviceV;` for moving from the device)
template <typename T, typename Y>
inline void
plmdDataToGPU (thrust::device_vector<T> &dvmem, Y &data, cudaStream_t stream) {
  plmdDataToGPU (dvmem, DataInterface (data), stream);
}

/// @brief copies data from the GPU, using thrust::device_vector as interface
///
/// if you are transferring data from/to a 'standard' vector, consider using
/// thrust::host_vector<T> and the thrust syntax (`deviceV=hostV;` for moving to
/// the device or `hostV=`deviceV;` for moving from the device)
template <typename T, typename Y>
inline void plmdDataFromGPU (thrust::device_vector<T> &dvmem, Y &data) {
  plmdDataFromGPU (dvmem, DataInterface (data));
}

/// @brief async version of plmdDataFromGPU
///
/// if you are transferring data from/to a 'standard' ( with a standar type,
/// like a int) vector, consider using thrust::host_vector<T> and the thrust
/// syntax (`deviceV=hostV;` for moving to the device or `hostV=`deviceV;` for
/// moving from the device)
template <typename T, typename Y>
inline void plmdDataFromGPU (thrust::device_vector<T> &dvmem,
                             Y &data,
                             cudaStream_t stream) {
  plmdDataFromGPU (dvmem, DataInterface (data), stream);
}

/// We'll find the ideal number of blocks using the Brent's theorem
size_t idealGroups (const size_t numberOfElements, const size_t runningThreads);

size_t threadsPerBlock (const unsigned N, const unsigned maxNumThreads);

/**********************************REDUCTIONS**********************************/
template <typename calculateFloat, int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void
reduction1DKernel (int num_valid,        // number if elements to be reduced
                   calculateFloat *d_in, // Tile of input
                   calculateFloat *d_out // Tile aggregate
                  ) {
  // Specialize BlockReduce type for our thread block
  using BlockReduceT = cub::BlockReduce<calculateFloat, BLOCK_THREADS>;
  // Shared memory
  __shared__ typename BlockReduceT::TempStorage temp_storage;
  const int data_id = threadIdx.x + blockIdx.x * blockDim.x;
  // Per-thread tile data
  calculateFloat data[ITEMS_PER_THREAD];
  cub::LoadDirectBlocked (data_id, d_in, data, num_valid, calculateFloat (0.0));
  // Compute sum
  calculateFloat aggregate = BlockReduceT (temp_storage).Sum (data);
  if (threadIdx.x == 0) {
    d_out[blockIdx.x] = aggregate;
  }
}

// the order of the template arguments is not standard: in this way T is deduced
// and the user can simply specify DATAPERTHREAD
template <unsigned DATAPERTHREAD, typename T, unsigned THREADS = 1024>
void doReduction1D_t (T *inputArray,
                      T *outputArray,
                      const unsigned int len,
                      const unsigned blocks,
                      const unsigned nthreads) {
  if constexpr (THREADS > 16) {
    // by using this "if constexpr" I do not need to add a specialized
    // declaration to end the loop
    if (nthreads == THREADS) {
      reduction1DKernel<T, THREADS, DATAPERTHREAD>
      <<<blocks, THREADS, THREADS * sizeof (T)>>> (
        len, inputArray, outputArray);
    } else {
      doReduction1D_t<DATAPERTHREAD, T, THREADS / 2> (
        inputArray, outputArray, len, blocks, nthreads);
    }
  } else {
    plumed_merror (
      "Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

template <unsigned DATAPERTHREAD, typename T, unsigned THREADS = 1024>
void doReduction1D_t (T *inputArray,
                      T *outputArray,
                      const unsigned int len,
                      const unsigned blocks,
                      const unsigned nthreads,
                      cudaStream_t stream) {
  if constexpr (THREADS > 16) {
    // by using this "if constexpr" I do not need to add a specialized
    // declaration to end the loop
    if (nthreads == THREADS) {
      reduction1DKernel<T, THREADS, DATAPERTHREAD>
      <<<blocks, THREADS, THREADS * sizeof (T), stream>>> (
        len, inputArray, outputArray);
    } else {
      doReduction1D_t<DATAPERTHREAD, T, THREADS / 2> (
        inputArray, outputArray, len, blocks, nthreads, stream);
    }
  } else {
    plumed_merror (
      "Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

template <unsigned DATAPERTHREAD, typename T>
void doReduction1D (T *inputArray,
                    T *outputArray,
                    const size_t len,
                    const unsigned blocks,
                    const size_t nthreads) {

  doReduction1D_t<DATAPERTHREAD> (
    inputArray, outputArray, len, blocks, nthreads);
}

template <unsigned DATAPERTHREAD, typename T>
void doReduction1D (T *inputArray,
                    T *outputArray,
                    const size_t len,
                    const unsigned blocks,
                    const size_t nthreads,
                    cudaStream_t stream) {

  doReduction1D_t<DATAPERTHREAD> (
    inputArray, outputArray, len, blocks, nthreads, stream);
}

template <typename calculateFloat, int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void
reductionNDKernel (int num_valid,        // number if elements to be reduced
                   calculateFloat *d_in, // Tile of input
                   calculateFloat *d_out // Tile aggregate
                  ) {
  // Specialize BlockReduce type for our thread block
  using BlockReduceT = cub::BlockReduce<calculateFloat, BLOCK_THREADS>;
  // Shared memory
  __shared__ typename BlockReduceT::TempStorage temp_storage;
  const int data_id = threadIdx.x + blockIdx.x * blockDim.x;
  calculateFloat data[ITEMS_PER_THREAD];
  cub::LoadDirectBlocked (data_id,
                          d_in + blockIdx.y * num_valid,
                          data,
                          num_valid,
                          calculateFloat (0.0));
  // Compute sum
  calculateFloat aggregate = BlockReduceT (temp_storage).Sum (data);
  if (threadIdx.x == 0) {
    d_out[blockIdx.x + blockIdx.y * gridDim.x] = aggregate;
  }
}

// the order of the template arguments is not standard: in this way T is deduced
// and the user can simply specify DATAPERTHREAD
template <unsigned DATAPERTHREAD, typename T, unsigned THREADS = 1024>
void doReductionND_t (T *inputArray,
                      T *outputArray,
                      const unsigned int len,
                      const dim3 blocks,
                      const unsigned nthreads) {
  if constexpr (THREADS > 16) {
    // by using this "if constexpr" I do not need to add a specialized
    // declaration to end the loop
    if (nthreads == THREADS) {
      reductionNDKernel<T, THREADS, DATAPERTHREAD>
      <<<blocks, THREADS, THREADS * sizeof (T)>>> (
        len, inputArray, outputArray);
    } else {
      doReductionND_t<DATAPERTHREAD, T, THREADS / 2> (
        inputArray, outputArray, len, blocks, nthreads);
    }
  } else {
    plumed_merror (
      "Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

template <unsigned DATAPERTHREAD, typename T, unsigned THREADS = 1024>
void doReductionND_t (T *inputArray,
                      T *outputArray,
                      const unsigned int len,
                      const dim3 blocks,
                      const unsigned nthreads,
                      cudaStream_t stream) {
  if constexpr (THREADS > 16) {
    // by using this "if constexpr" I do not need to add a specialized
    // declaration to end the loop
    if (nthreads == THREADS) {
      reductionNDKernel<T, THREADS, DATAPERTHREAD>
      <<<blocks, THREADS, THREADS * sizeof (T), stream>>> (
        len, inputArray, outputArray);
    } else {
      doReductionND_t<DATAPERTHREAD, T, THREADS / 2> (
        inputArray, outputArray, len, blocks, nthreads, stream);
    }
  } else {
    plumed_merror (
      "Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

template <unsigned DATAPERTHREAD, typename T>
void doReductionND (T *inputArray,
                    T *outputArray,
                    const unsigned int len,
                    const dim3 blocks,
                    const unsigned nthreads) {

  doReductionND_t<DATAPERTHREAD> (
    inputArray, outputArray, len, blocks, nthreads);
}

template <unsigned DATAPERTHREAD, typename T>
void doReductionND (T *inputArray,
                    T *outputArray,
                    const unsigned int len,
                    const dim3 blocks,
                    const unsigned nthreads,
                    cudaStream_t stream) {

  doReductionND_t<DATAPERTHREAD> (
    inputArray, outputArray, len, blocks, nthreads, stream);
}
} // namespace CUDAHELPERS
#endif //__PLUMED_cuda_helpers_cuh
