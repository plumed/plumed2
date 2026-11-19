/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 Reactive Voronoi contributors
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
// FP64 streamed assignments with optional buffered center lists.
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <new>

namespace {
constexpr int denseThreads=128;
constexpr int neighborThreads=32;

bool checkCuda(cudaError_t status) {
  if(status==cudaSuccess) {
    return true;
  }
  std::fprintf(stderr,"ReactiveVoronoi CUDA: %s\n",cudaGetErrorString(status));
  return false;
}

/// Preserve the device context selected by the MD engine around every call.
class DeviceScope {
  int previous_=-1;
  bool valid_=false;
public:
  explicit DeviceScope(int device) {
    int current=-1;
    if(!checkCuda(cudaGetDevice(&current))) {
      return;
    }
    if(current!=device) {
      if(!checkCuda(cudaSetDevice(device))) {
        return;
      }
      previous_=current;
    }
    valid_=true;
  }
  ~DeviceScope() {
    if(previous_>=0) {
      cudaSetDevice(previous_);
    }
  }
  DeviceScope(const DeviceScope&)=delete;
  DeviceScope& operator=(const DeviceScope&)=delete;
  bool valid() const {
    return valid_;
  }
};

struct State {
  int numCenters=0;
  int numAssigned=0;
  int device=-1;
  double* positions=nullptr;
  double* occupancy=nullptr;
  double* maximum=nullptr;
  double* normalization=nullptr;
  double* sensitivity=nullptr;
  double* derivatives=nullptr;
  double* virial=nullptr;
  int* invalid=nullptr;
  int* offsets=nullptr;
  int* candidateCenters=nullptr;
  int candidateCapacity=0;
  double cutoff=0;
  bool useNeighbors=false;
  double box[3] {};
  double kappa=0;
  bool pbc=false;
  cudaStream_t stream=nullptr;
};

template<int threads>
__device__ double reduce(double value, bool maximum, double* scratch) {
  const int thread=threadIdx.x;
  if constexpr(threads==neighborThreads) {
    for(int offset=threads/2; offset; offset/=2) {
      const double other=__shfl_down_sync(0xffffffff,value,offset);
      value=maximum ? fmax(value,other) : value+other;
    }
    return __shfl_sync(0xffffffff,value,0);
  }
  scratch[thread]=value;
  __syncthreads();
  for(int offset=threads/2; offset; offset/=2) {
    if(thread<offset) {
      scratch[thread]=maximum ?
                      fmax(scratch[thread],scratch[thread+offset]) :
                      scratch[thread]+scratch[thread+offset];
    }
    __syncthreads();
  }
  const double result=scratch[0];
  __syncthreads();
  return result;
}

__device__ double distance(const State state, int center, int assigned,
                           double* displacement) {
  double squared=0;
  for(int axis=0; axis<3; ++axis) {
    displacement[axis]=state.positions[3*(state.numCenters+assigned)+axis]-
                       state.positions[3*center+axis];
    if(state.pbc) {
      displacement[axis]-=nearbyint(displacement[axis]/state.box[axis])*state.box[axis];
    }
    squared+=displacement[axis]*displacement[axis];
  }
  return sqrt(squared);
}

template<int threads, bool useNeighbors>
__global__ void forwardKernel(State state) {
  __shared__ double scratch[threads];
  const int assigned=blockIdx.x;
  const int thread=threadIdx.x;
  const int begin=useNeighbors ? state.offsets[assigned] : 0;
  const int end=useNeighbors ? state.offsets[assigned+1] : state.numCenters;
  double maximum=-INFINITY;
  for(int entry=begin+thread; entry<end; entry+=threads) {
    const int center=useNeighbors ? state.candidateCenters[entry] : entry;
    double displacement[3];
    const double length=distance(state,center,assigned,displacement);
    if(!isfinite(length) || length<=2.220446049250313e-16) {
      atomicExch(state.invalid,1);
    }
    if(!useNeighbors || length<=state.cutoff) {
      const double score=-state.kappa*length;
      if(!isfinite(score)) {
        atomicExch(state.invalid,1);
      }
      maximum=fmax(maximum,score);
    }
  }
  maximum=reduce<threads>(maximum,true,scratch);
  if(useNeighbors && maximum==-INFINITY) {
    if(thread==0) {
      atomicCAS(state.invalid,0,2);
    }
    return;
  }
  double normalization=0;
  for(int entry=begin+thread; entry<end; entry+=threads) {
    const int center=useNeighbors ? state.candidateCenters[entry] : entry;
    double displacement[3];
    const double length=distance(state,center,assigned,displacement);
    if(!useNeighbors || length<=state.cutoff) {
      normalization+=exp(-state.kappa*length-maximum);
    }
  }
  normalization=reduce<threads>(normalization,false,scratch);
  if(thread==0) {
    state.maximum[assigned]=maximum;
    state.normalization[assigned]=normalization;
  }
  for(int entry=begin+thread; entry<end; entry+=threads) {
    const int center=useNeighbors ? state.candidateCenters[entry] : entry;
    double displacement[3];
    const double length=distance(state,center,assigned,displacement);
    if(useNeighbors && length>state.cutoff) {
      continue;
    }
    const double weight=exp(-state.kappa*length-maximum)/normalization;
    atomicAdd(state.occupancy+center,weight);
  }
}

template<int threads, bool useNeighbors>
__global__ void backwardKernel(State state) {
  __shared__ double scratch[threads];
  const int assigned=blockIdx.x;
  const int thread=threadIdx.x;
  const int begin=useNeighbors ? state.offsets[assigned] : 0;
  const int end=useNeighbors ? state.offsets[assigned+1] : state.numCenters;
  const double maximum=state.maximum[assigned];
  const double normalization=state.normalization[assigned];
  double mean=0;
  for(int entry=begin+thread; entry<end; entry+=threads) {
    const int center=useNeighbors ? state.candidateCenters[entry] : entry;
    double displacement[3];
    const double length=distance(state,center,assigned,displacement);
    if(useNeighbors && length>state.cutoff) {
      continue;
    }
    const double weight=exp(-state.kappa*length-maximum)/normalization;
    mean+=weight*state.sensitivity[center];
  }
  mean=reduce<threads>(mean,false,scratch);
  double assignedDerivative[3]= {0,0,0};
  double virial[9]= {};
  for(int entry=begin+thread; entry<end; entry+=threads) {
    const int center=useNeighbors ? state.candidateCenters[entry] : entry;
    double displacement[3];
    const double length=distance(state,center,assigned,displacement);
    if(useNeighbors && length>state.cutoff) {
      continue;
    }
    const double weight=exp(-state.kappa*length-maximum)/normalization;
    const double radial=-state.kappa*weight*(state.sensitivity[center]-mean)/length;
    for(int axis=0; axis<3; ++axis) {
      const double gradient=radial*displacement[axis];
      atomicAdd(state.derivatives+3*center+axis,-gradient);
      assignedDerivative[axis]+=gradient;
      for(int component=0; component<3; ++component) {
        virial[3*axis+component]-=gradient*displacement[component];
      }
    }
  }
  for(int axis=0; axis<3; ++axis) {
    const double value=reduce<threads>(assignedDerivative[axis],false,scratch);
    if(thread==0) {
      state.derivatives[3*(state.numCenters+assigned)+axis]=value;
    }
  }
  for(int component=0; component<9; ++component) {
    const double value=reduce<threads>(virial[component],false,scratch);
    if(thread==0) {
      atomicAdd(state.virial+component,value);
    }
  }
}

void release(State* state) {
  if(!state) {
    return;
  }
  DeviceScope scope(state->device);
  if(scope.valid()) {
    cudaFree(state->positions);
    cudaFree(state->occupancy);
    cudaFree(state->maximum);
    cudaFree(state->normalization);
    cudaFree(state->sensitivity);
    cudaFree(state->derivatives);
    cudaFree(state->virial);
    cudaFree(state->invalid);
    cudaFree(state->offsets);
    cudaFree(state->candidateCenters);
    if(state->stream) {
      cudaStreamDestroy(state->stream);
    }
  }
  delete state;
}
}

extern "C" int reactive_voronoi_device_count_v1() {
  int count=0;
  if(cudaGetDeviceCount(&count)!=cudaSuccess) {
    return 0;
  }
  return count;
}

extern "C" void* reactive_voronoi_create_v2(int numCenters, int numAssigned,
    int mdDevice) {
  const int deviceCount=reactive_voronoi_device_count_v1();
  if(numCenters<=0 || numAssigned<=0 || deviceCount==0 || mdDevice>=deviceCount) {
    return nullptr;
  }
  int device=mdDevice;
  if(device<0 && !checkCuda(cudaGetDevice(&device))) {
    return nullptr;
  }
  DeviceScope scope(device);
  if(!scope.valid()) {
    return nullptr;
  }
  State* state=new(std::nothrow) State;
  if(!state) {
    return nullptr;
  }
  state->numCenters=numCenters;
  state->numAssigned=numAssigned;
  state->device=device;
  const size_t numAtoms=static_cast<size_t>(numCenters)+numAssigned;
  if(!checkCuda(cudaStreamCreateWithFlags(&state->stream,cudaStreamNonBlocking)) ||
      !checkCuda(cudaMalloc(&state->positions,3*numAtoms*sizeof(double))) ||
      !checkCuda(cudaMalloc(&state->occupancy,numCenters*sizeof(double))) ||
      !checkCuda(cudaMalloc(&state->maximum,numAssigned*sizeof(double))) ||
      !checkCuda(cudaMalloc(&state->normalization,numAssigned*sizeof(double))) ||
      !checkCuda(cudaMalloc(&state->sensitivity,numCenters*sizeof(double))) ||
      !checkCuda(cudaMalloc(&state->derivatives,3*numAtoms*sizeof(double))) ||
      !checkCuda(cudaMalloc(&state->virial,9*sizeof(double))) ||
      !checkCuda(cudaMalloc(&state->invalid,sizeof(int)))) {
    release(state);
    return nullptr;
  }
  return state;
}

extern "C" void* reactive_voronoi_create_v1(int numCenters, int numAssigned) {
  return reactive_voronoi_create_v2(numCenters,numAssigned,-1);
}

extern "C" int reactive_voronoi_device_v1(void* context) {
  return context ? static_cast<State*>(context)->device : -1;
}

extern "C" void reactive_voronoi_destroy_v1(void* context) {
  release(static_cast<State*>(context));
}

namespace {
int runForward(void* context, const double* positions, const double* box, int pbc,
               double kappa, double* occupancy, bool useNeighbors, double cutoff,
               const int* offsets, const int* centers, int count) {
  if(!context || !positions || !box || !occupancy || !std::isfinite(kappa) || kappa<=0) {
    return 1;
  }
  State& state=*static_cast<State*>(context);
  DeviceScope scope(state.device);
  if(!scope.valid()) {
    return 1;
  }
  if(useNeighbors) {
    if(offsets) {
      if(count<0 || (count>0 && !centers)) {
        return 1;
      }
      if(!state.offsets && !checkCuda(cudaMalloc(&state.offsets,
                                      (static_cast<size_t>(state.numAssigned)+1)*sizeof(int)))) {
        return 1;
      }
      if(count>state.candidateCapacity) {
        int* replacement=nullptr;
        if(!checkCuda(cudaMalloc(&replacement,static_cast<size_t>(count)*sizeof(int)))) {
          return 1;
        }
        cudaFree(state.candidateCenters);
        state.candidateCenters=replacement;
        state.candidateCapacity=count;
      }
      if(!checkCuda(cudaMemcpyAsync(state.offsets,offsets,
                                    (static_cast<size_t>(state.numAssigned)+1)*sizeof(int),
                                    cudaMemcpyHostToDevice,state.stream)) ||
          (count>0 && !checkCuda(cudaMemcpyAsync(state.candidateCenters,centers,
                                 static_cast<size_t>(count)*sizeof(int),cudaMemcpyHostToDevice,state.stream)))) {
        return 1;
      }
    } else if(!state.offsets) {
      return 1;
    }
    state.cutoff=cutoff;
  }
  state.useNeighbors=useNeighbors;
  state.pbc=pbc;
  state.kappa=kappa;
  for(int axis=0; axis<3; ++axis) {
    state.box[axis]=box[axis];
    if(pbc && (!std::isfinite(box[axis]) || box[axis]<=0)) {
      return 1;
    }
  }
  const size_t numAtoms=static_cast<size_t>(state.numCenters)+state.numAssigned;
  if(!checkCuda(cudaMemcpyAsync(state.positions,positions,3*numAtoms*sizeof(double),
                                cudaMemcpyHostToDevice,state.stream)) ||
      !checkCuda(cudaMemsetAsync(state.occupancy,0,state.numCenters*sizeof(double),state.stream)) ||
      !checkCuda(cudaMemsetAsync(state.invalid,0,sizeof(int),state.stream))) {
    return 1;
  }
  if(useNeighbors) {
    forwardKernel<neighborThreads,true><<<state.numAssigned,neighborThreads,0,state.stream>>>(state);
  } else {
    forwardKernel<denseThreads,false><<<state.numAssigned,denseThreads,0,state.stream>>>(state);
  }
  int invalid=0;
  if(!checkCuda(cudaGetLastError()) ||
      !checkCuda(cudaMemcpyAsync(occupancy,state.occupancy,state.numCenters*sizeof(double),
                                 cudaMemcpyDeviceToHost,state.stream)) ||
      !checkCuda(cudaMemcpyAsync(&invalid,state.invalid,sizeof(int),
                                 cudaMemcpyDeviceToHost,state.stream)) ||
      !checkCuda(cudaStreamSynchronize(state.stream))) {
    return 1;
  }
  return invalid;
}

}

extern "C" int reactive_voronoi_forward_v1(void* context, const double* positions,
    const double* box, int pbc, double kappa, double* occupancy) {
  return runForward(context,positions,box,pbc,kappa,occupancy,false,0,nullptr,nullptr,0);
}

// A null offsets pointer reuses the already uploaded candidate list.
extern "C" int reactive_voronoi_forward_nlist_v1(void* context, const double* positions,
    const double* box, int pbc, double kappa, double* occupancy, double cutoff,
    const int* offsets, const int* centers, int count) {
  return runForward(context,positions,box,pbc,kappa,occupancy,true,cutoff,offsets,centers,count);
}

extern "C" int reactive_voronoi_backward_v1(void* context, const double* sensitivity,
    double* derivatives, double* virial) {
  if(!context || !sensitivity || !derivatives || !virial) {
    return 1;
  }
  State& state=*static_cast<State*>(context);
  DeviceScope scope(state.device);
  if(!scope.valid()) {
    return 1;
  }
  const size_t numAtoms=static_cast<size_t>(state.numCenters)+state.numAssigned;
  if(!checkCuda(cudaMemcpyAsync(state.sensitivity,sensitivity,state.numCenters*sizeof(double),
                                cudaMemcpyHostToDevice,state.stream)) ||
      !checkCuda(cudaMemsetAsync(state.derivatives,0,3*numAtoms*sizeof(double),state.stream)) ||
      !checkCuda(cudaMemsetAsync(state.virial,0,9*sizeof(double),state.stream))) {
    return 1;
  }
  if(state.useNeighbors) {
    backwardKernel<neighborThreads,true><<<state.numAssigned,neighborThreads,0,state.stream>>>(state);
  } else {
    backwardKernel<denseThreads,false><<<state.numAssigned,denseThreads,0,state.stream>>>(state);
  }
  if(!checkCuda(cudaGetLastError()) ||
      !checkCuda(cudaMemcpyAsync(derivatives,state.derivatives,3*numAtoms*sizeof(double),
                                 cudaMemcpyDeviceToHost,state.stream)) ||
      !checkCuda(cudaMemcpyAsync(virial,state.virial,9*sizeof(double),
                                 cudaMemcpyDeviceToHost,state.stream)) ||
      !checkCuda(cudaStreamSynchronize(state.stream))) {
    return 1;
  }
  return 0;
}
