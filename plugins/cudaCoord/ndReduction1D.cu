#include "ndReduction.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cudaHelpers.cuh"

//std::min is constexpr in c++14
//#include <algorithm>

#include <iomanip>
#include <iostream>
#include <vector>
#include <numeric>
#define vdbg(...) std::cerr << std::setw(4) << __LINE__ <<":" << std::setw(20)<< #__VA_ARGS__ << " " << (__VA_ARGS__) <<'\n'
//#define vdbg(...)

// * Grids map to GPUs
// * Blocks map to the MultiProcessors (MP)
// * Threads map to Stream Processors (SP)
// * Warps are groups of (32) threads that execute simultaneously


//There are a LOTS of unrolled loop down here,
//see this https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf to undestand why

//https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#dim3
//7.3.2. dim3
//This type is an integer vector type based on uint3 that is used to specify dimensions.

namespace PLMD {
namespace CUDAHELPERS {

template <unsigned numThreads, typename T>
__device__ void warpReduce1D(volatile T* sdata, unsigned int place) {
  if(numThreads >= 64) { //compile time
    sdata[place] += sdata[place + 32];
  }
  if(numThreads >= 32) { //compile time
    sdata[place] += sdata[place + 16];
  }
  if(numThreads >= 16) { //compile time
    sdata[place] += sdata[place + 8];
  }
  if(numThreads >= 8) { //compile time
    sdata[place] += sdata[place + 4];
  }
  if(numThreads >= 4) { //compile time
    sdata[place] += sdata[place + 2];
  }
  if(numThreads >= 2) { //compile time
    sdata[place] += sdata[place + 1];
  }
}

template <unsigned numThreads, typename T>
__global__ void reduction1D(T *g_idata, T *g_odata, const unsigned int len) {
  //extern __shared__ T sdata[numThreads];
  auto sdata = shared_memory_proxy<T>();
  const unsigned int place = threadIdx.x;
  // each thread loads one element from global to shared mem
  unsigned int i = numThreads*blockIdx.x*2 + place;
  const unsigned int gridSize = numThreads*gridDim.x*2;
  sdata[place] = T(0);
  //I think this may slow down the loop, but this does not force the user to have
  //an input that is multiple of the threads, padded with zeros
  while (i+numThreads < len) {
    sdata[place] += g_idata[i] + g_idata[i+numThreads];
    i+=gridSize;
  }
  while (i < len) {
    sdata[place] += g_idata[i];
    i+=gridSize;
  }

  __syncthreads();
  // do reduction in shared memory

  if (numThreads >= 512) {//compile time
    if (threadIdx.x  < 256) {
      sdata[place] += sdata[place + 256];
    } __syncthreads();
  }
  if (numThreads >= 256) {//compile time
    if (threadIdx.x  < 128) {
      sdata[place] += sdata[place + 128];
    } __syncthreads();
  }
  if (numThreads >= 128) {//compile time
    if (threadIdx. x < 64) {
      sdata[place] += sdata[place + 64];
    } __syncthreads();
  }
  //Instructions are SIMD synchronous within a warp
  //so no need for __syncthreads(), in the last iterations
  if (threadIdx.x < mymin(32u,numThreads/2)) {
    warpReduce1D<numThreads>(sdata, place);
  }
  // write result for this block to global mem
  if (threadIdx.x == 0) {
    g_odata[blockIdx.x] = sdata[0];
  }
}

//after c++14 the template activation will be shorter to write:
//template<typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>

///finds the nearest upper multiple of the given reference (wit non increments)
template<typename T,
         typename std::enable_if<std::is_integral<T>::value, bool>::type = true>
inline T nearestUpperMultipleTo(T number, T reference) {
  return ((number-1)|(reference-1))+1;
}

///We'll find the ideal number of blocks using the Brent's theorem
size_t getIdealGroups(size_t numberOfElements, size_t runningThreads) {
  //nearest upper multiple to the numberof threads
  const size_t nnToGPU=nearestUpperMultipleTo(numberOfElements,runningThreads);
  ///Brent’s theorem says each thread should sum O(log n) elements
  //const size_t elementsPerThread=log(nnToGPU);
  const size_t expectedTotalThreads = ceil(nnToGPU/log(nnToGPU));
  //hence the blocks should have this size:
  const unsigned ngroups = nearestUpperMultipleTo(expectedTotalThreads,runningThreads)/runningThreads;
  return  ngroups;
}

template <typename T>
void callReduction1D (T *g_idata, T *g_odata, const unsigned int len, const unsigned blocks, const unsigned nthreads) {
  switch (nthreads) {
  case 512:
    reduction1D<512,T><<<blocks,512,512*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 256:
    reduction1D<256,T><<<blocks,256,256*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 128:
    reduction1D<128,T><<<blocks,128,128*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 64:
    reduction1D<64, T><<<blocks,64,64*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 32:
    reduction1D<32, T><<<blocks,32,32*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  }
}

size_t decideThreadsPerBlock(unsigned N, unsigned maxNumThreads=512) {
  //this seeks the minimum number of threads to use a sigle block (and end the recursion)
  size_t dim=32;
  for (dim=32; dim<512; dim<<=1) {
    if (maxNumThreads < dim) {
      dim >>=1;
      break;
    }
    if( N < dim) {
      break;
    }
  }
  return dim;
}

double reduceScalar(double* cudaScalarAddress, unsigned N, unsigned maxNumThreads) {
//we'll proceed to call recursively callreduction1D until N==1:
  double *reduceOut = cudaScalarAddress;
  double *reduceIn;
  vdbg("In");
  while(N>1) {
    size_t runningThreads = decideThreadsPerBlock(N,maxNumThreads);
    reduceIn = reduceOut;
    reduceOut = nullptr;
    vdbg(N);
    vdbg(cudaScalarAddress);
    vdbg(reduceIn);
    vdbg(reduceOut);
    vdbg(runningThreads);
    auto ngroups=getIdealGroups(N, runningThreads);
    vdbg(ngroups);
    cudaMalloc(&reduceOut,ngroups  * sizeof(double));
    vdbg(reduceOut);
    callReduction1D (reduceIn, reduceOut, N, ngroups, runningThreads);
    if (reduceIn != cudaScalarAddress) {
      vdbg("Free reduceIn");
      cudaFree(reduceIn);
    }
    N=ngroups;
  }
  double toret;
  cudaMemcpy(&toret, reduceOut, sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(reduceOut);
  vdbg(toret);
  return toret;
}
} //namespace CUDAHELPERS
} //namespace PLMD
/**todo:
 * @compiletime request all the possible threads 1 2 4 8 16 32 64 128 256 512
 * create a function that cases the threadsnum
 * pass an already initializated cudavector to be reduced with the needed dimensions
 *
*/

