#include "ndReduction.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iomanip>
#include <iostream>
#include <vector>
#include <numeric>

//#define vdbg(...) std::cerr << std::setw(4) << __LINE__ <<":" << std::setw(20)<< #__VA_ARGS__ << " " << (__VA_ARGS__) <<'\n'
#define vdbg(...)

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
__global__ void reductionND(const T *g_idata, T *g_odata, const unsigned int len) {
  //playing with this 
  //https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
  auto sdata = shared_memory_proxy<T>();
  const unsigned int coord = blockIdx.y;
  const unsigned int place = threadIdx.x;
  // each thread loads one element from global to shared memory
  const unsigned int diplacement = blockIdx.y*len;
  unsigned int i = (numThreads*2)*blockIdx.x + place + diplacement;
  const unsigned int gridSize = (numThreads*2)*gridDim.x;
  const unsigned int trgt=len+diplacement;

  sdata[threadIdx.x] = T(0);
  while (i+numThreads < trgt) {
    sdata[threadIdx.x] += g_idata[i] + g_idata[i+numThreads];
    i+=gridSize;
  }
  while (i < trgt) {
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
    g_odata[blockIdx.x+blockIdx.y*gridDim.x] = sdata[0];    
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

template <typename T>
void callReductionND (T *g_idata, T *g_odata, const unsigned int len, const dim3 blocks, const unsigned nthreads){
    switch (nthreads) {
  case 512:
    reductionND<512,T><<<blocks,512,512*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 256:
    reductionND<256,T><<<blocks,256,256*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 128:
    reductionND<128,T><<<blocks,128,128*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 64:
    reductionND<64, T><<<blocks,64,64*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  case 32:
    reductionND<32, T><<<blocks,32,32*sizeof(T)>>>(g_idata,g_odata, len);
    break;
  default:
    plumed_merror("Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
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
    cudaFree(reduceOut);
    cudaMalloc(&reduceOut,ngroups  * sizeof(double));
    callReduction1D (reduceIn, reduceOut, N, ngroups, runningThreads);
    if (reduceIn != cudaScalarAddress){
      cudaFree(reduceIn);
    }
    N=ngroups;
  }
  double toret;
  cudaMemcpy(&toret, reduceOut, sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(reduceOut);
  return toret;
}

double reduceScalar(memoryHolder<double>& cudaScalarAddress,
 memoryHolder<double>& memoryHelper,
  unsigned N, unsigned maxNumThreads){
//we'll proceed to call recursively callreduction1D until N==1:
  memoryHolder<double>* reduceIn= &memoryHelper;
  memoryHolder<double>* reduceOut =&cudaScalarAddress;
  while(N>1){
    size_t runningThreads = decideThreadsPerBlock(N,maxNumThreads);
    std::swap(reduceIn,reduceOut);
    auto ngroups=getIdealGroups(N, runningThreads);
    reduceOut->resize(ngroups);
    callReduction1D (reduceIn->getPointer(), reduceOut->getPointer(), N, ngroups, runningThreads);
    N=ngroups;
  }
  double toret;
  reduceOut->copyFromCuda(&toret);
  return toret;
}

std::vector<Vector> reduceNVectors(double* cudaNVectorAddress, unsigned N, unsigned nat, unsigned maxNumThreads){
  double *reduceOut = cudaNVectorAddress;
  double *reduceIn;
  auto dim = nat*3;
  vdbg("InNVectors");
  while(N>1){
    size_t runningThreads = decideThreadsPerBlock(N,maxNumThreads);
    reduceIn = reduceOut;
    reduceOut = nullptr;
    vdbg(N);
    vdbg(cudaNVectorAddress);
    vdbg(reduceIn);
    vdbg(reduceOut);
    vdbg(runningThreads);
    dim3 ngroups(getIdealGroups(N, runningThreads),dim);
    vdbg(ngroups.x);
    vdbg(ngroups.y);
    cudaFree(reduceOut);
    cudaMalloc(&reduceOut,ngroups.y* ngroups.x  * sizeof(double));
    vdbg(reduceOut);

    callReductionND (reduceIn, reduceOut, N, ngroups, runningThreads);

    if (reduceIn != cudaNVectorAddress){
      vdbg("Free reduceIn");
      cudaFree(reduceIn);
    }
    N=ngroups.x;
  }
  std::vector<Vector> toret(nat);
  cudaMemcpy(&toret[0][0], reduceOut, 3*nat*sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(reduceOut);
  vdbg(toret[0]);
  return toret;
}

//THIS DOES NOT KEEP THE DATA SAFE
std::vector<Vector> reduceNVectors(memoryHolder<double>& cudaNVectorAddress,
 memoryHolder<double>& memoryHelper, 
unsigned N, unsigned nat, unsigned maxNumThreads){
  memoryHolder<double>* reduceIn= &memoryHelper;
  memoryHolder<double>* reduceOut =&cudaNVectorAddress;
  
  auto dim = nat*3;
  vdbg("InNVectors");
  while(N>1){
    size_t runningThreads = decideThreadsPerBlock(N,maxNumThreads);
    std::swap(reduceIn,reduceOut);
    dim3 ngroups(getIdealGroups(N, runningThreads),dim);
    reduceOut->resize(ngroups.y* ngroups.x);
    

    callReductionND (reduceIn->getPointer(), reduceOut->getPointer(), N, ngroups, runningThreads);

    N=ngroups.x;
  }
  std::vector<Vector> toret(nat);
  reduceOut->copyFromCuda(&toret[0][0]);
  //cudaMemcpy(&toret[0][0], reduceOut, 3*nat*sizeof(double), cudaMemcpyDeviceToHost);
  
  vdbg(toret[0]);
  return toret;
}

Vector reduceVector(double* cudaVectorAddress, unsigned N, unsigned maxNumThreads){
///@TODO:This is not tested as now
//we'll proceed to call recursively callreduction1D until N==1:
  double *reduceOut = cudaVectorAddress;
  double *reduceIn;
  while(N>1){
    size_t runningThreads = decideThreadsPerBlock(N,maxNumThreads);
    reduceIn = reduceOut;
    reduceOut = nullptr;
    dim3 ngroups(getIdealGroups(N, runningThreads),3);
    cudaFree(reduceOut);
    cudaMalloc(&reduceOut,ngroups.y* ngroups.x  * sizeof(double));
    
    callReductionND (reduceIn, reduceOut, N, ngroups, runningThreads);
        
    if (reduceIn != cudaVectorAddress){
      cudaFree(reduceIn);
    }
    N=ngroups.x;
  }
  Vector toret;
  cudaMemcpy(&toret[0], reduceOut, 3*sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(reduceOut);
  return toret;
}

//#define vdbg(...) std::cerr << std::setw(4) << __LINE__ <<":" << std::setw(20)<< #__VA_ARGS__ << " " << (__VA_ARGS__) <<'\n'
Tensor reduceTensor(memoryHolder<double>&  cudaTensorAddress, 
memoryHolder<double>& memoryHelper, unsigned N, unsigned maxNumThreads){
//we'll proceed to call recursively callreduction1D until N==1:
  memoryHolder<double>* reduceIn= &memoryHelper;
  memoryHolder<double>* reduceOut =&cudaTensorAddress;
  while(N>1){
    size_t runningThreads = decideThreadsPerBlock(N,maxNumThreads);
    std::swap(reduceIn,reduceOut);
    dim3 ngroups(getIdealGroups(N, runningThreads),9);
    reduceOut->resize(ngroups.y* ngroups.x);
  
    callReductionND (reduceIn->getPointer(), reduceOut->getPointer(),
    N, ngroups, runningThreads);
        
    N=ngroups.x;
  }
  Tensor toret;
  reduceOut->copyFromCuda(&toret[0][0]);
  return toret;
}

DVS::DVS(unsigned nat): deriv(nat){}

template <typename T>
void callReduction1D (T *g_idata, T *g_odata, const unsigned int len,
 const unsigned blocks, const unsigned nthreads,cudaStream_t& stream){
  switch (nthreads) {
  case 512:
    reduction1D<512,T><<<blocks,512,512*sizeof(T),stream>>>(g_idata,g_odata, len);
    break;
  case 256:
    reduction1D<256,T><<<blocks,256,256*sizeof(T),stream>>>(g_idata,g_odata, len);
    break;
  case 128:
    reduction1D<128,T><<<blocks,128,128*sizeof(T),stream>>>(g_idata,g_odata, len);
    break;
  case 64:
    reduction1D<64, T><<<blocks,64,64*sizeof(T),stream>>>(g_idata,g_odata, len);
    break;
  case 32:
    reduction1D<32, T><<<blocks,32,32*sizeof(T),stream>>>(g_idata,g_odata, len);
    break;
  default:
    plumed_merror("Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

template <typename T>
void callReductionND (T *g_idata, T *g_odata, const unsigned int len,
 const dim3 blocks, const unsigned nthreads,cudaStream_t& stream){
    switch (nthreads) {
  case 512:
    reductionND<512,T><<<blocks,512,512*sizeof(T),stream>>>(g_idata,g_odata, len);
    break;
  case 256:
    reductionND<256,T><<<blocks,256,256*sizeof(T),stream>>>(g_idata,g_odata, len);
    break;
  case 128:
    reductionND<128,T><<<blocks,128,128*sizeof(T),stream>>>(g_idata,g_odata, len);
    break;
  case 64:
    reductionND<64, T><<<blocks,64,64*sizeof(T),stream>>>(g_idata,g_odata, len);
    break;
  case 32:
    reductionND<32, T><<<blocks,32,32*sizeof(T),stream>>>(g_idata,g_odata, len);
    break;
  default:
    plumed_merror("Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}


DVS reduceDVS(memoryHolder<double>& cudaV,
memoryHolder<double>& cudaT,
memoryHolder<double>& cudaS,
 memoryHolder<double>& memoryHelperV, 
 memoryHolder<double>& memoryHelperT, 
 memoryHolder<double>& memoryHelperS, 
unsigned N, unsigned nat, unsigned maxNumThreads){
  memoryHolder<double>* reduceSIn= &memoryHelperS;
  memoryHolder<double>* reduceSOut =&cudaS;
  memoryHolder<double>* reduceTIn= &memoryHelperT;
  memoryHolder<double>* reduceTOut =&cudaT;
  memoryHolder<double>* reduceVIn= &memoryHelperV;
  memoryHolder<double>* reduceVOut =&cudaV;
  
  auto dim = nat*3;
  vdbg("InNVectors");
  cudaStream_t streamV;
  cudaStream_t streamT;
  cudaStream_t streamS;
  cudaStreamCreate(&streamV);
  cudaStreamCreate(&streamT);
  cudaStreamCreate(&streamS);
  while(N>1){
    vdbg(N);
    size_t runningThreads = decideThreadsPerBlock(N,maxNumThreads);
    std::swap(reduceTIn,reduceTOut);
    std::swap(reduceSIn,reduceSOut);
    std::swap(reduceVIn,reduceVOut);
    unsigned ngroupsS=getIdealGroups(N, runningThreads);
    dim3 ngroupsV(ngroupsS,dim);
    dim3 ngroupsT(ngroupsS,9);
    
    reduceVOut->resize(ngroupsV.y* ngroupsV.x);
    reduceTOut->resize(ngroupsT.y* ngroupsT.x);
    reduceSOut->resize(ngroupsS);
    vdbg("Calls");
    callReductionND (reduceVIn->getPointer(), reduceVOut->getPointer(), N, ngroupsV, runningThreads,streamV);
    vdbg("Calls");
    callReductionND (reduceTIn->getPointer(), reduceTOut->getPointer(), N, ngroupsT, runningThreads,streamT);
    vdbg("Calls");
    callReduction1D (reduceSIn->getPointer(), reduceSOut->getPointer(), N, ngroupsS, runningThreads,streamS);
    
    cudaDeviceSynchronize();
    vdbg("End");
    N=ngroupsS;
  }
  DVS toret(nat);
  reduceVOut->copyFromCuda(&toret.deriv[0][0]);
  reduceTOut->copyFromCuda(&toret.virial[0][0]);
  reduceSOut->copyFromCuda(&toret.scalar);
  cudaStreamDestroy(streamV);
  cudaStreamDestroy(streamT);
  cudaStreamDestroy(streamS);
  return toret;
}


} //namespace CUDAHELPERS
} //namespace PLMD

