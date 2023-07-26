#include "ndReduction.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cudaHelpers.cuh"

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
__device__ void warpReduce(volatile T* sdata, unsigned int place){
    if(numThreads >= 64){//compile time
      sdata[place] += sdata[place + 32];
    }
    if(numThreads >= 32){//compile time
      sdata[place] += sdata[place + 16];
    }
    if(numThreads >= 16){//compile time
      sdata[place] += sdata[place + 8];
    }
    if(numThreads >= 8){//compile time
      sdata[place] += sdata[place + 4];
    }
    if(numThreads >= 4){//compile time
      sdata[place] += sdata[place + 2];
    }
    if(numThreads >= 2){//compile time
      sdata[place] += sdata[place + 1];
    }
}

template <unsigned numThreads, typename T>
__global__ void reductionND(T *g_idata, T *g_odata, const unsigned int len) {
  //playing with this 
  //https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
  auto sdata = shared_memory_proxy<T>();
  const unsigned int coord = blockIdx.y;
  const unsigned int dim = gridDim.y;
  const unsigned int place = dim*threadIdx.x+coord;
  const unsigned int totalLen = dim*len;
  // each thread loads one element from global to shared mem
  //const unsigned int tid = dim*threadIdx.x+coord;
  unsigned int i = (numThreads*2)*blockIdx.x*dim + place;
  const unsigned int gridSize = (numThreads*2)*gridDim.x*dim;
  sdata[threadIdx.x] = T(0);
  while (i+numThreads*dim < totalLen) {
    sdata[threadIdx.x] += g_idata[i] + g_idata[i+numThreads*dim];    
    i+=gridSize;
  }
  while (i < totalLen) {
    sdata[threadIdx.x] += g_idata[i];
     i+=gridSize;
  }

  __syncthreads();
  // do reduction in shared memory
  
  if (numThreads >= 512) {//compile time
    if (threadIdx.x  < 256) {
       sdata[threadIdx.x] += sdata[threadIdx.x + 256]; } __syncthreads(); 
    }
  if (numThreads >= 256) {//compile time
    if (threadIdx.x  < 128) {
       sdata[threadIdx.x] += sdata[threadIdx.x + 128]; } __syncthreads(); 
    }
  if (numThreads >= 128) {//compile time
    if (threadIdx. x < 64) { 
      sdata[threadIdx.x] += sdata[threadIdx.x + 64]; } __syncthreads();
    }
  //Instructions are SIMD synchronous within a warp
  //so no need for __syncthreads(), in the last iterations
  if (threadIdx.x < mymin(32u,numThreads/2)) {
    warpReduce<numThreads>(sdata, threadIdx.x);
  }
  // write result for this block to global memory
  if (threadIdx.x == 0){
    g_odata[dim*blockIdx.x+coord] = sdata[0];
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
       sdata[place] += sdata[place + 256]; } __syncthreads(); 
       }
  if (numThreads >= 256) {//compile time
    if (threadIdx.x  < 128) {
       sdata[place] += sdata[place + 128]; } __syncthreads(); 
       }
  if (numThreads >= 128) {//compile time
    if (threadIdx. x < 64) { 
      sdata[place] += sdata[place + 64]; } __syncthreads();
       }
  //Instructions are SIMD synchronous within a warp
  //so no need for __syncthreads(), in the last iterations
  if (threadIdx.x < mymin(32u,numThreads/2)) {
    warpReduce<numThreads>(sdata, place);
  }
  // write result for this block to global mem
  if (threadIdx.x == 0){
    g_odata[blockIdx.x] = sdata[0];
  }
}

//after c++14 the template activation will be shorter to write:
//template<typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>

///finds the nearest upper multiple of the given reference (wit non increments)
template<typename T, 
typename std::enable_if<std::is_integral<T>::value, bool>::type = true>
  inline T nearestUpperMultipleTo(T number, T reference){
    return ((number-1)|(reference-1))+1;
}

///We'll find the ideal number of blocks using the Brent's theorem
size_t getIdealGroups(size_t numberOfElements, size_t runningThreads){
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
void callReduction1D (T *g_idata, T *g_odata, const unsigned int len, const unsigned blocks, const unsigned nthreads){
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
  default:
    plumed_merror("Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

size_t decideThreadsPerBlock(unsigned N, unsigned maxNumThreads=512){
  //this seeks the minimum number of threads to use a sigle block (and end the recursion)
  size_t dim=32;
  for (dim=32;dim<512;dim<<=1){
    if (maxNumThreads < dim) {
      dim >>=1;
      break;
    }
    if( N < dim){
      break;
    }
  }
  return dim;
}

double reduceScalar(double* cudaScalarAddress, unsigned N, unsigned maxNumThreads){
//we'll proceed to call recursively callreduction1D until N==1:
  double *reduceOut = cudaScalarAddress;
  double *reduceIn;
  while(N>1){
    size_t runningThreads = decideThreadsPerBlock(N,maxNumThreads);
    reduceIn = reduceOut;
    reduceOut = nullptr;
    auto ngroups=getIdealGroups(N, runningThreads);
    cudaMalloc(&reduceOut,ngroups  * sizeof(double));
    callReduction1D (reduceIn, reduceOut, N, ngroups, runningThreads);
    /*
    {
    std::vector<double> coordsToSUM(N);
    cudaMemcpy(coordsToSUM.data(), reduceIn, N*sizeof(double), cudaMemcpyDeviceToHost);
    double ncoordIn=std::accumulate(coordsToSUM.begin(),coordsToSUM.end(),0.0);
    vdbg(ncoordIn);
    }
    {
    std::vector<double> coordsToSUM(ngroups);
    cudaMemcpy(coordsToSUM.data(), reduceOut, ngroups*sizeof(double), cudaMemcpyDeviceToHost);
    double ncoordOut=std::accumulate(coordsToSUM.begin(),coordsToSUM.end(),0.0);
    vdbg(ncoordOut);
    }
    */
    if (reduceIn != cudaScalarAddress){
      cudaFree(reduceIn);
    }
    N=ngroups;
  }
  double toret;
  cudaMemcpy(&toret, reduceOut, sizeof(double), cudaMemcpyDeviceToHost);
  return toret;
}

template <typename T>
void callReductionND (T *g_idata, T *g_odata, const unsigned int len, const dim3 blocks, const unsigned nthreads){
  const unsigned N = blocks.y;
  switch (nthreads) {
  case 512:
    reductionND<512,T><<<blocks,512,N*512*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 256:
    reductionND<256,T><<<blocks,256,N*256*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 128:
    reductionND<128,T><<<blocks,128,N*128*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 64:
    reductionND<64, T><<<blocks,64,N*64*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 32:
    reductionND<32, T><<<blocks,32,N*32*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  default:
    plumed_merror("Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

Tensor reduceTensor(double* cudaTensorAddress, unsigned N, unsigned maxNumThreads){
//we'll proceed to call recursively callreduction1D until N==1:
  double *reduceOut = cudaTensorAddress;
  double *reduceIn;
  vdbg("InTensor");
  while(N>1){
    size_t runningThreads = decideThreadsPerBlock(N,maxNumThreads);
    reduceIn = reduceOut;
    reduceOut = nullptr;
    vdbg(N);
    vdbg(cudaTensorAddress);
    vdbg(reduceIn);
    vdbg(reduceOut);
    vdbg(runningThreads);
    dim3 ngroups(getIdealGroups(N, runningThreads),9);
    vdbg(ngroups.x);
    vdbg(ngroups.y);
    cudaMalloc(&reduceOut,ngroups.y* ngroups.x  * sizeof(double));
    vdbg(reduceOut);    
        {
    std::vector<double> coordsToSUM(ngroups.y* N);
    cudaMemcpy(coordsToSUM.data(), reduceIn, ngroups.y* N *sizeof(double), cudaMemcpyDeviceToHost);
    for (auto j=0u; j< ngroups.y; ++j){
      double t=0.0;
    for (auto i=0u; i< N; ++i){
      t+= coordsToSUM[i*ngroups.y+j];
    }
    std::cerr << __LINE__ << " ti= " << t <<"\n";
    }
    
    }
    callReductionND (reduceIn, reduceOut, N, ngroups, runningThreads);
    {
    std::vector<double> coordsToSUM(ngroups.y* ngroups.x);
    cudaMemcpy(coordsToSUM.data(), reduceOut, ngroups.y* ngroups.x *sizeof(double), cudaMemcpyDeviceToHost);
    for (auto j=0u; j< ngroups.y; ++j){
      double t=0.0;
    for (auto i=0u; i< ngroups.x; ++i){
      t+= coordsToSUM[i*ngroups.y+j];
    }
    std::cerr << __LINE__ << " to= " << t <<"\n";
    }
    
    }
    
    if (reduceIn != cudaTensorAddress){
      vdbg("Free reduceIn");
      cudaFree(reduceIn);
    }
    N=ngroups.x;
  }
  Tensor toret;
  cudaMemcpy(&toret[0][0], reduceOut, 9*sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(reduceOut);
  vdbg(toret);
  return toret;
}
} //namespace CUDAHELPERS
} //namespace PLMD

