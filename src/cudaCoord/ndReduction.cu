#include "ndReduction.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iomanip>
#include <iostream>
#include <vector>
#include <numeric>

//#define vdbg(...) std::cerr << std::setw(4) << __LINE__ <<":" << std::setw(20)<< #__VA_ARGS__ << " " << (__VA_ARGS__) <<'\n'
#define vdbg(...)

// Some help for me:
// * Grids map to GPUs
// * Blocks map to the MultiProcessors (MP)
// * Threads map to Stream Processors (SP)
// * Warps are groups of (32) threads that execute simultaneously


//There are a LOTS of unrolled loop down here,
//see this https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
// to undestand why

//https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#dim3
//7.3.2. dim3
//This type is an integer vector type based on uint3 that is used to specify dimensions. 

namespace PLMD {
namespace CUDAHELPERS {
template <unsigned numThreads, typename T>
__global__ void reductionND(const T *inputArray, T *outputArray, const unsigned int len) {
  //playing with this 
  //https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
  auto sdata = shared_memory_proxy<T>();
  //const unsigned int coord = blockIdx.y;
  const unsigned int place = threadIdx.x;
  // each thread loads one element from global to shared memory
  const unsigned int diplacement = blockIdx.y*len;
  unsigned int i = (numThreads*2)*blockIdx.x + place + diplacement;
  const unsigned int gridSize = (numThreads*2)*gridDim.x;
  //the first element is in blockIdx.y*len, the last element to sum in (blockIdx.y+1)*len-1
  const unsigned int trgt=len+diplacement;

  sdata[threadIdx.x] = T(0);
  while (i+numThreads < trgt) {
    sdata[threadIdx.x] += inputArray[i] + inputArray[i+numThreads];
    i+=gridSize;
  }
  while (i < trgt) {
    sdata[threadIdx.x] += inputArray[i];
     i+=gridSize;
  }

  __syncthreads();
  // do reduction in shared memory
  reductor<numThreads>(sdata,outputArray,blockIdx.x+blockIdx.y*gridDim.x);
}

template <unsigned numThreads, typename T>
__global__ void reduction1D(T *inputArray, T *outputArray, const unsigned int len) {
  auto sdata = shared_memory_proxy<T>();
  const unsigned int place = threadIdx.x;
  // each thread sums some elements from global to shared memory
  unsigned int i = (2*numThreads)*blockIdx.x + place;
  const unsigned int gridSize = (2*numThreads)*gridDim.x;
  sdata[place] = T(0);
  //The double while is for preventig wrong memory
  while ( i + numThreads < len) {
    sdata[place] += inputArray[i] + inputArray[i+numThreads];
    i+=gridSize;
  }
  while (i < len) {
    sdata[place] += inputArray[i];
    i+=gridSize;
  }
    
  __syncthreads();
  // do reduction in shared memory
  reductor<numThreads>(sdata,outputArray,blockIdx.x);
}

template <unsigned numThreads, typename T>
__global__ void reductionDerivatives(T *g_idata, T *g_odata,unsigned* nnlist, const unsigned int len) {
  //we saved the dd in an array x0 x1 x2..,xn-1,y0 y1 y2..,yn-1,z0 z1 z2..,zn-1
  //virialOut[ii*3+jj]-=d[ii]*d[jj]*dfunc;
  //printf("CUDA:****\n");
  auto sdata = shared_memory_proxy<T>();
  //const unsigned int xyz = blockIdx.z;
  const unsigned int place = threadIdx.x;
  const unsigned int atomId = blockIdx.y;
  // each thread loads one element from global to shared mem
  unsigned int nn = (2*numThreads) * blockIdx.x + place;
  const unsigned int i = /*nn +*/ blockIdx.z * len;
  const unsigned int dfunc = /*nn +*/ 3 * len;
  const unsigned int gridSize = (2*numThreads) * gridDim.x;

  sdata[place] = T(0);
  //I think this may slow down the loop, but this does not force the user to have
  //an input that is multiple of the threads, padded with zeros
  while (nn + numThreads < len) {
    
    if(atomId == nnlist[nn*2]){
      sdata[place] -= g_idata[nn+i]*g_idata[nn+dfunc];
    //} else if(atomId==nnlist[nn*2+1]){
      sdata[place] += g_idata[nn+i]*g_idata[nn+dfunc];
    //}
    
    //if(atomId == nnlist[(nn+numThreads)*2]){
      sdata[place] -= g_idata[nn+i+numThreads]*g_idata[nn+dfunc+numThreads];
    //} else if (atomId == nnlist[(nn+numThreads)*2+1]){
      sdata[place]  += g_idata[nn+i+numThreads]*g_idata[nn+dfunc+numThreads];
    //}
    //   g_idata[i]*g_idata[dfunc]
    // + g_idata[i+numThreads]*g_idata[dfunc+numThreads];
    //i+=gridSize;
    
    nn+=gridSize;
    //dfunc+=gridSize;
  }
  while (nn < len) {
    if(atomId == nnlist[nn*2]){
      sdata[place] -= g_idata[nn+i]*g_idata[nn+dfunc];
    } else if(atomId==nnlist[nn*2+1]){
      sdata[place] += g_idata[nn+i]*g_idata[nn+dfunc];
    }
     
    //i+=gridSize;
    nn+=gridSize;
    //dfunc+=gridSize;
  }
    
  __syncthreads();
  // do reduction in shared memory
  reductor<numThreads>(sdata,g_odata,
    blockIdx.x //position
    + 3 * gridDim.x * blockIdx.y // dispacement by atom number 
    + gridDim.x * blockIdx.z //coordinate (x,y,z) displacement
  );
  //   if(atomId==(gridDim.y-1)
  // //  &&
  // //    blockIdx.x==0
  //  && 
  //  threadIdx.x==0
  // ) {
  //     printf("CUDA:****  %i %s at=%i: %f\n",blockIdx.x,
  //       ((blockIdx.z==0)?"x":((blockIdx.z==1)?"y":"z")),atomId,sdata[0]);
  //     //printf("CUDA:****  x:%i y:%i z:%i\n",gridDim.x,gridDim.y,gridDim.z);
  //  }
}

//after c++14 the template activation will be shorter to write:
//template<typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>

///finds the nearest upper multiple of the given reference
template<typename T, 
typename std::enable_if<std::is_integral<T>::value, bool>::type = true>
  inline T nearestUpperMultipleTo(T number, T reference){
    return ((number-1)|(reference-1))+1;
}

///We'll find the ideal number of blocks using the Brent's theorem
size_t idealGroups(size_t numberOfElements, size_t runningThreads){
    //nearest upper multiple to the numberof threads
    const size_t nnToGPU=nearestUpperMultipleTo(numberOfElements,runningThreads);
    ///Brent’s theorem says each thread should sum O(log n) elements
    //const size_t elementsPerThread=log(nnToGPU);
    const size_t expectedTotalThreads = ceil(nnToGPU/log(nnToGPU));
    //hence the blocks should have this size:
    const unsigned ngroups = nearestUpperMultipleTo(expectedTotalThreads,runningThreads)/runningThreads;
    return  ngroups;
}


size_t threadsPerBlock(unsigned N, unsigned maxNumThreads){
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
void callReduction1D (
  T *inputArray,
  T *outputArray,
  const unsigned int len,
  const unsigned blocks,
  const unsigned nthreads){
  switch (nthreads) {
  case 512:
    reduction1D<512,T><<<blocks,512,512*sizeof(T)>>>(inputArray,outputArray, len);
    break;
  case 256:
    reduction1D<256,T><<<blocks,256,256*sizeof(T)>>>(inputArray,outputArray, len);
    break;
  case 128:
    reduction1D<128,T><<<blocks,128,128*sizeof(T)>>>(inputArray,outputArray, len);
    break;
  case 64:
    reduction1D<64, T><<<blocks,64,64*sizeof(T)>>>(inputArray,outputArray, len);
    break;
  case 32:
    reduction1D<32, T><<<blocks,32,32*sizeof(T)>>>(inputArray,outputArray, len);
    break;
  default:
    plumed_merror("Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

template <typename T>
void callReductionND (T *inputArray,
  T *outputArray,
  const unsigned int len,
  const dim3 blocks,
  const unsigned nthreads){
    switch (nthreads) {
  case 512:
    reductionND<512,T><<<blocks,512,512*sizeof(T)>>>(inputArray,outputArray, len);
    break;
  case 256:
    reductionND<256,T><<<blocks,256,256*sizeof(T)>>>(inputArray,outputArray, len);
    break;
  case 128:
    reductionND<128,T><<<blocks,128,128*sizeof(T)>>>(inputArray,outputArray, len);
    break;
  case 64:
    reductionND<64, T><<<blocks,64,64*sizeof(T)>>>(inputArray,outputArray, len);
    break;
  case 32:
    reductionND<32, T><<<blocks,32,32*sizeof(T)>>>(inputArray,outputArray, len);
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
    size_t runningThreads = threadsPerBlock(N,maxNumThreads);
    reduceIn = reduceOut;
    reduceOut = nullptr;
    auto ngroups=idealGroups(N, runningThreads);
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
    size_t runningThreads = threadsPerBlock(N,maxNumThreads);
    std::swap(reduceIn,reduceOut);
    auto ngroups=idealGroups(N, runningThreads);
    reduceOut->resize(ngroups);
    callReduction1D (reduceIn->pointer(), reduceOut->pointer(), N, ngroups, runningThreads);
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
  while(N>1){
    size_t runningThreads = threadsPerBlock(N,maxNumThreads);
    reduceIn = reduceOut;
    reduceOut = nullptr;
    dim3 ngroups(idealGroups(N, runningThreads),dim);
    cudaFree(reduceOut);
    cudaMalloc(&reduceOut,ngroups.y* ngroups.x  * sizeof(double));
    callReductionND (reduceIn, reduceOut, N, ngroups, runningThreads);

    if (reduceIn != cudaNVectorAddress){
      cudaFree(reduceIn);
    }
    N=ngroups.x;
  }
  std::vector<Vector> toret(nat);
  cudaMemcpy(&toret[0][0], reduceOut, 3*nat*sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(reduceOut);
  return toret;
}

//THIS DOES NOT KEEP THE DATA SAFE
std::vector<Vector> reduceNVectors(
  memoryHolder<double>& cudaNVectorAddress,
  memoryHolder<double>& memoryHelper, 
  unsigned N,
  unsigned nat,
  unsigned maxNumThreads){
  memoryHolder<double>* reduceIn= &memoryHelper;
  memoryHolder<double>* reduceOut =&cudaNVectorAddress;
  
  auto dim = nat*3;
  while(N>1){
    size_t runningThreads = threadsPerBlock(N,maxNumThreads);
    std::swap(reduceIn,reduceOut);
    dim3 ngroups(idealGroups(N, runningThreads),dim);
    reduceOut->resize(ngroups.y* ngroups.x);
    callReductionND (
      reduceIn->pointer(),
      reduceOut->pointer(),
      N,
      ngroups,
      runningThreads);
    N=ngroups.x;
  }
  std::vector<Vector> toret(nat);
  reduceOut->copyFromCuda(&toret[0][0]);
  //cudaMemcpy(&toret[0][0], reduceOut, 3*nat*sizeof(double), cudaMemcpyDeviceToHost);
  return toret;
}

Vector reduceVector(double* cudaVectorAddress, unsigned N, unsigned maxNumThreads){
///@TODO:This is not tested as now
//we'll proceed to call recursively callreduction1D until N==1:
  double *reduceOut = cudaVectorAddress;
  double *reduceIn;
  while(N>1){
    size_t runningThreads = threadsPerBlock(N,maxNumThreads);
    reduceIn = reduceOut;
    reduceOut = nullptr;
    dim3 ngroups(idealGroups(N, runningThreads),3);
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

Tensor reduceTensor(memoryHolder<double>&  cudaTensorAddress, 
memoryHolder<double>& memoryHelper, unsigned N, unsigned maxNumThreads){
//we'll proceed to call recursively callreduction1D until N==1:
  memoryHolder<double>* reduceIn= &memoryHelper;
  memoryHolder<double>* reduceOut =&cudaTensorAddress;
  while(N>1){
    size_t runningThreads = threadsPerBlock(N,maxNumThreads);
    std::swap(reduceIn,reduceOut);
    dim3 ngroups(idealGroups(N, runningThreads),9);
    reduceOut->resize(ngroups.y* ngroups.x);
  
    callReductionND (reduceIn->pointer(), reduceOut->pointer(),
    N, ngroups, runningThreads);
        
    N=ngroups.x;
  }
  Tensor toret;
  reduceOut->copyFromCuda(&toret[0][0]);
  return toret;
}

template <typename T>
void callReduction1D (T *inputArray, T *outputArray, const unsigned int len,
 const unsigned blocks, const unsigned nthreads,cudaStream_t stream=0){
  switch (nthreads) {
  case 512:
    reduction1D<512,T><<<blocks,512,512*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 256:
    reduction1D<256,T><<<blocks,256,256*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 128:
    reduction1D<128,T><<<blocks,128,128*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 64:
    reduction1D<64, T><<<blocks,64,64*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 32:
    reduction1D<32, T><<<blocks,32,32*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  default:
    plumed_merror("Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

template <typename T>
void callReductionND (T *inputArray, T *outputArray, const unsigned int len,
 const dim3 blocks, const unsigned nthreads,cudaStream_t stream=0){
    switch (nthreads) {
  case 512:
    reductionND<512,T><<<blocks,512,512*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 256:
    reductionND<256,T><<<blocks,256,256*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 128:
    reductionND<128,T><<<blocks,128,128*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 64:
    reductionND<64, T><<<blocks,64,64*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  case 32:
    reductionND<32, T><<<blocks,32,32*sizeof(T),stream>>>(inputArray,outputArray, len);
    break;
  default:
    plumed_merror("Reduction can be called only with 512, 256, 128, 64 or 32 threads.");
  }
}

void doReduction1D (double *inputArray,
 double *outputArray,
 const unsigned int len,
 const unsigned blocks,
 const unsigned nthreads,
 cudaStream_t stream){
  callReduction1D (inputArray, outputArray, len, blocks, nthreads, stream);
 }

void doReductionND (double *inputArray,
 double *outputArray,
 const unsigned int len,
 const dim3 blocks,
 const unsigned nthreads,
 cudaStream_t stream){
  callReductionND (inputArray, outputArray, len, blocks, nthreads, stream);
 }

 void doReduction1D (float *inputArray,
 float *outputArray,
 const unsigned int len,
 const unsigned blocks,
 const unsigned nthreads,
 cudaStream_t stream){
  callReduction1D (inputArray, outputArray, len, blocks, nthreads, stream);
 }

void doReductionND (float *inputArray,
 float *outputArray,
 const unsigned int len,
 const dim3 blocks,
 const unsigned nthreads,
 cudaStream_t stream){
  callReductionND (inputArray, outputArray, len, blocks, nthreads, stream);
 }
} // namespace CUDAHELPERS
} // namespace PLMD
