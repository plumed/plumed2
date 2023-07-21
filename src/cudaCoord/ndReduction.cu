#include <vector>

#include <random>
#include <iostream>


#include <cassert>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define vdbg(...) std::cerr << __LINE__ <<":" << #__VA_ARGS__ << " " << (__VA_ARGS__) <<'\n'
//#define vdbg(...)

using rndengine=std::mt19937;

template <unsigned numThreads, typename T>
__device__ void warpReduceND(volatile T* sdata, const unsigned int place, const unsigned int dim){
    if(numThreads >= 64){//compile time
      sdata[place] += sdata[place + dim*32];
      }
    if(numThreads >= 32){//compile time
      sdata[place] += sdata[place + dim*16];
      }
    if(numThreads >= 16){//compile time
      sdata[place] += sdata[place + dim*8];
      }
    if(numThreads >= 8){//compile time
      sdata[place] += sdata[place + dim*4];
      }
    if(numThreads >= 4){//compile time
      sdata[place] += sdata[place + dim*2];
      }
    if(numThreads >= 2){//compile time
      sdata[place] += sdata[place + dim];
      }
}

template <unsigned numThreads, typename T>
__global__ void reductionND(T *g_idata, T *g_odata, const unsigned int len) {
  //playing with this 
  //https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
  extern __shared__ T sdata[];
  const unsigned int coord = blockIdx.y;
  const unsigned int dim = gridDim.y;
  const unsigned int place = dim*threadIdx.x+coord;
  // each thread loads one element from global to shared mem
  //const unsigned int tid = dim*threadIdx.x+coord;
  unsigned int i = (numThreads*2)*blockIdx.x*dim + place;
  const unsigned int gridSize = (numThreads*2)*gridDim.x*dim;
  sdata[place] = T(0);
  while (i < len) {
    sdata[place] += g_idata[i] + g_idata[i+numThreads*dim];
    i+=gridSize;
  }
    
  __syncthreads();
  // do reduction in shared memory
  
  if (numThreads >= 512) {//compile time
    if (threadIdx.x  < 256) {
       sdata[place] += sdata[place + 256*dim]; } __syncthreads(); 
       }
  if (numThreads >= 256) {//compile time
    if (threadIdx.x  < 128) {
       sdata[place] += sdata[place + 128*dim]; } __syncthreads(); 
       }
  if (numThreads >= 128) {//compile time
    if (threadIdx. x < 64) { 
      sdata[place] += sdata[place + 64*dim]; } __syncthreads();
       }
  //Instructions are SIMD synchronous within a warp
  //so no need for __syncthreads(), in the last iterations
  if (threadIdx.x < 32) {
    warpReduceND<numThreads>(sdata, place,dim);
  }
  // write result for this block to global mem
  if (threadIdx.x == 0){
    //printf("thread [%i],%i: %i\n",blockIdx.x, threadIdx.y,sdata[threadIdx.y] );
    g_odata[dim*blockIdx.x+coord] = sdata[coord];
  }
}

template <unsigned numThreads, typename T>
__device__ void warpReduce1D(volatile T* sdata, unsigned int place){
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
__global__ void reduction1D(T *g_idata, T *g_odata, const unsigned int len) {
  extern __shared__ T sdata[numThreads];
  const unsigned int place = threadIdx.x;
  // each thread loads one element from global to shared mem
  unsigned int i = numThreads*blockIdx.x*2 + place;
  const unsigned int gridSize = numThreads*gridDim.x*2;
  sdata[place] = T(0);
  while (i < len) {
    sdata[place] += g_idata[i] + g_idata[i+numThreads];
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
  if (threadIdx.x < 32) {
    warpReduce1D<numThreads>(sdata, place);
  }
  // write result for this block to global mem
  if (threadIdx.x == 0){
    g_odata[blockIdx.x] = sdata[0];
  }
}

//c++14 ->template<typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>

///finds the nearest upper multiple of the given reference
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

template <unsigned numThreads, typename T>
std::vector<T> reductionCuda1D(const std::vector<T>& in){
    const auto size = in.size();
    const unsigned ngroups = getIdealGroups(size,numThreads);
  
    T * cudaCoords;
    //no need to alloc nnToGPU things, (see getIdealGroups)
    // because we are passing size as limiting item
    cudaMalloc(&cudaCoords,  size * sizeof(T));
    //cudaMemset(cudaCoords, 0,nnToGPU * sizeof(T));
    cudaMemcpy(cudaCoords, in.data(), size * sizeof(T), cudaMemcpyHostToDevice);
    T * cudaRet;
    cudaMalloc(&cudaRet,ngroups  * sizeof(T));
    
    reduction1D<numThreads><<<ngroups, numThreads>>>(cudaCoords,cudaRet,size);
      
    std::vector<T> toret(ngroups,0);
    cudaMemcpy(toret.data(), cudaRet, ngroups * sizeof(T), cudaMemcpyDeviceToHost);
    cudaFree(cudaRet);
    cudaFree(cudaCoords);
    if (ngroups == 1){
      return toret;
    }
    return reductionCuda1D<numThreads>(toret);
}

template <unsigned numThreads, typename T>
std::vector<T> reductionCuda(const std::vector<T>& in, unsigned N){
  if (N==1){
    return reductionCuda1D<numThreads>(in);
  }

  //TODO: limit the memory consumptions and the thread number
    const auto size = in.size();
    if (size%N != 0){
        throw "Please give a ND Vector";
    }
    const auto rsize=size/N;
    const unsigned ngroups = getIdealGroups(rsize,numThreads);
    
    const auto retSize = ngroups *N;

    T * cudaCoords;
    cudaMalloc(&cudaCoords, N * size * sizeof(T));
    cudaMemcpy(cudaCoords, in.data(), size * sizeof(T), cudaMemcpyHostToDevice);
    T * cudaRet;
    cudaMalloc(&cudaRet,ngroups * N * sizeof(T));
    cudaMemset(cudaRet, 0, N * sizeof(T));
    dim3 groups2D(ngroups,N);
    
    reductionND<numThreads><<<groups2D, numThreads, N*numThreads*sizeof(T)>>>(cudaCoords,cudaRet,size);
    
    std::vector<T> toret(retSize,0);
    cudaMemcpy(toret.data(), cudaRet, retSize * sizeof(T), cudaMemcpyDeviceToHost);
    cudaFree(cudaRet);
    cudaFree(cudaCoords);
    if (retSize == N){
      return toret;
    }
    
    return reductionCuda<numThreads>(toret,  N);
      
}

void tmp(){
  
std::vector<int> in(0);
reductionCuda<1>(in, 1);
reductionCuda1D<1>(in);
}

/**todo:
 * @compiletaime request all the possible threads 1 2 4 8 16 32 64 128 256 512
 * create a function that cases the threadsnum
 * pass an already initializated cudavector to be reduced with the needed dimensions
 * 
*/