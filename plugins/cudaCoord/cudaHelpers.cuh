#ifndef __PLUMED_cuda_helpers_cuh
#define __PLUMED_cuda_helpers_cuh

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <vector>
namespace PLMD {
namespace CUDAHELPERS {

template <class T> __device__ constexpr const T &mymin(const T &a, const T &b) {
  return (b < a) ? b : a;
}

template <typename T> __device__ T *shared_memory_proxy() {
  // do we need an __align__() here?
  extern __shared__ unsigned char memory[];
  return reinterpret_cast<T *>(memory);
}

// onlymoveable helper
template <typename T> class memoryHolder {
  /*
  struct info {
    T* pointer_{nullptr};
    //the dimension in memory
    unsigned dim_{0};
    //the used part of the array
    unsigned usedim_{0};
    //mutex
    //number of implementations
  };
  */
  T *pointer_{nullptr};
  // the dimension in memory
  unsigned dim_{0};
  // the used part of the array
  unsigned usedim_{0};

public:
  memoryHolder() {}
  memoryHolder(memoryHolder &other) = delete;
  memoryHolder(memoryHolder &&other)
    : pointer_(other.pointer_), dim_(other.dim_), usedim_(other.usedim_) {
    other.pointer_ = nullptr;
  }
  memoryHolder(const unsigned newDim) : dim_(newDim), usedim_(newDim) {
    cudaMalloc(&pointer_, dim_ * sizeof(T));
  }
  memoryHolder &operator=(memoryHolder &&other) {
    pointer_ = other.pointer_;
    dim_ = other.dim_;
    usedim_ = other.usedim_;
    other.pointer_ = nullptr;
    other.dim_ = 0;
    other.usedim_ = 0;
    return *this;
  }

  ~memoryHolder() { cudaFree(pointer_); }
  T *pointer() { return pointer_; }
  unsigned size() const { return usedim_; }
  unsigned reserved() const { return dim_; }
  void /*CUresult*/ resize(const unsigned newDim, bool keepData = false) {

    if (newDim > dim_) {
      T *new_ptr = nullptr;
      cudaMalloc(&new_ptr, newDim * sizeof(T));
      if (keepData && pointer_) {
        // copy the old data into the new pointer
        cudaMemcpy(new_ptr, pointer_, usedim_ * sizeof(T),
                   cudaMemcpyDeviceToDevice);
      }
      cudaFree(pointer_);
      pointer_ = new_ptr;
      dim_ = newDim;
    }
    usedim_ = newDim;
  }
  void /*CUresult*/ reserve(const unsigned newDim, bool keepData = false) {
    if (newDim > dim_) {
      T *new_ptr = nullptr;
      cudaMalloc(&new_ptr, newDim * sizeof(T));
      if (keepData && pointer_) {
        // copy the old data into the new pointer
        cudaMemcpy(new_ptr, pointer_, usedim_ * sizeof(T));
      }
      cudaFree(pointer_);
      pointer_ = new_ptr;
      dim_ = newDim;
    }
  }

  void /*CUresult*/ copyToCuda(T *cpupointer) {
    cudaMemcpy(pointer_, cpupointer, usedim_ * sizeof(T),
               cudaMemcpyHostToDevice);
  }

  void /*CUresult*/ copyToCuda(T *cpupointer, cudaStream_t stream) {
    cudaMemcpyAsync(pointer_, cpupointer, usedim_ * sizeof(T),
                    cudaMemcpyHostToDevice, stream);
  }

  void /*CUresult*/ copyFromCuda(T *cpupointer) {
    cudaMemcpy(cpupointer, pointer_, usedim_ * sizeof(T),
               cudaMemcpyDeviceToHost);
  }

  void /*CUresult*/ copyFromCuda(T *cpupointer, cudaStream_t stream) {
    cudaMemcpyAsync(cpupointer, pointer_, usedim_ * sizeof(T),
                    cudaMemcpyDeviceToHost, stream);
  }

  //*conversions*/
  // using vector to avoid memory leaks
  template <typename Y> void /*CUresult*/ copyToCuda(Y *cpupointer) {
    std::vector<T> tempMemory(usedim_);
    for (auto i = 0u; i < usedim_; ++i)
      tempMemory[i] = cpupointer[i];

    cudaMemcpy(pointer_, tempMemory.data(), usedim_ * sizeof(T),
               cudaMemcpyHostToDevice);
  }

  template <typename Y>
  void /*CUresult*/ copyToCuda(Y *cpupointer, cudaStream_t stream) {
    std::vector<T> tempMemory(usedim_);
    for (auto i = 0u; i < usedim_; ++i)
      tempMemory[i] = cpupointer[i];

    cudaMemcpyAsync(pointer_, tempMemory.data(), usedim_ * sizeof(T),
                    cudaMemcpyHostToDevice, stream);
  }

  template <typename Y> void /*CUresult*/ copyFromCuda(Y *cpupointer) {
    std::vector<T> tempMemory(usedim_);
    cudaMemcpy(tempMemory.data(), pointer_, usedim_ * sizeof(T),
               cudaMemcpyDeviceToHost);
    for (auto i = 0u; i < usedim_; ++i)
      cpupointer[i] = tempMemory[i];
  }

  template <typename Y>
  void /*CUresult*/ copyFromCuda(Y *cpupointer, cudaStream_t stream) {
    std::vector<T> tempMemory(usedim_);
    cudaMemcpyAsync(tempMemory.data(), pointer_, usedim_ * sizeof(T),
                    cudaMemcpyDeviceToHost, stream);
    for (auto i = 0u; i < usedim_; ++i)
      cpupointer[i] = tempMemory[i];
  }
  /**/

  void swap(memoryHolder &other) {
    std::swap(this->pointer_, other.pointer_);
    std::swap(this->dim_, other.dim_);
    std::swap(this->usedim_, other.usedim_);
  }
};

template <typename T> void swap(memoryHolder<T> &a, memoryHolder<T> &b) {
  // a swap with no no copies :)
  memoryHolder<T> t = std::move(a);
  a = std::move(b);
  b = std::move(t);
}

// template <unsigned numThreads, typename T>
// __device__ void warpReduce(volatile T* sdata, const unsigned int place){
//   //Instructions are SIMD synchronous within a warp
//   //so no need for __syncthreads(), in the last iterations
//     if(numThreads >= 64){//compile time
//       sdata[place] += sdata[place + 32];
//     }
//     if(numThreads >= 32){//compile time
//       sdata[place] += sdata[place + 16];
//     }
//     if(numThreads >= 16){//compile time
//       sdata[place] += sdata[place + 8];
//     }
//     if(numThreads >= 8){//compile time
//       sdata[place] += sdata[place + 4];
//     }
//     if(numThreads >= 4){//compile time
//       sdata[place] += sdata[place + 2];
//     }
//     if(numThreads >= 2){//compile time
//       sdata[place] += sdata[place + 1];
//     }
// }

template <unsigned numThreads, typename T>
__device__ void reductor(volatile T *sdata, T *outputArray,
                         const unsigned int where) {
  const unsigned int tid = threadIdx.x;
  if (numThreads >= 1024) { // compile time
    if (tid < 512)
      sdata[tid] += sdata[tid + 512];
    __syncthreads();
  }
  if (numThreads >= 512) { // compile time
    if (tid < 256)
      sdata[tid] += sdata[tid + 256];
    __syncthreads();
  }
  if (numThreads >= 256) { // compile time
    if (tid < 128)
      sdata[tid] += sdata[tid + 128];
    __syncthreads();
  }
  if (numThreads >= 128) { // compile time
    if (tid < 64)
      sdata[tid] += sdata[tid + 64];
    __syncthreads();
  }
  if (tid < mymin(32u, numThreads / 2)) {
    // warpReduce<numThreads>(sdata, tid);
    if (numThreads >= 64) { // compile time
      sdata[tid] += sdata[tid + 32];
    }
    if (numThreads >= 32) { // compile time
      sdata[tid] += sdata[tid + 16];
    }
    if (numThreads >= 16) { // compile time
      sdata[tid] += sdata[tid + 8];
    }
    if (numThreads >= 8) { // compile time
      sdata[tid] += sdata[tid + 4];
    }
    if (numThreads >= 4) { // compile time
      sdata[tid] += sdata[tid + 2];
    }
    if (numThreads >= 2) { // compile time
      sdata[tid] += sdata[tid + 1];
    }
  }
  // write result for this block to global memory
  if (tid == 0) {
    outputArray[where] = sdata[0];
  }
}

} // namespace CUDAHELPERS
} // namespace PLMD
#endif //__PLUMED_cuda_helpers_cuh
