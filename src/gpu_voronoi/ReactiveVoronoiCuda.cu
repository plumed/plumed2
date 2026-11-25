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
// FP64 soft assignment of the reactive Voronoi actions on the GPU (CUDA or HIP),
// with complete assignments or with the candidate list of NLIST.
#include "colvar/ReactiveVoronoiBackend.h"
#include "core/GpuDevice.h"
#include "tools/Pbc.h"
#include "tools/Tools.h"
#include "gpu_cuda_hip/cudaHelpers.cuh"
#include <algorithm>
#include <climits>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace PLMD {
namespace gpu_voronoi {
namespace {

using CUDAHELPERS::cudaCheck;
using CUDAHELPERS::deviceBuffer;
using CUDAHELPERS::pinnedBuffer;
using CUDAHELPERS::warpMask;

constexpr int denseThreads=128;
constexpr int neighborThreads=32;

/// What the kernels see: the device buffers and the parameters of the last forward pass
struct Data {
  int numCenters=0;
  int numAssigned=0;
  const double* positions=nullptr;
  double* occupancy=nullptr;
  double* maximum=nullptr;
  double* normalization=nullptr;
  const double* sensitivity=nullptr;
  double* derivatives=nullptr;
  double* virial=nullptr;
  unsigned long long* invalid=nullptr;
  const int* offsets=nullptr;
  const int* candidateCenters=nullptr;
  double cutoff=0;
  double box[3] {};
  double inverseBox[3] {};
  double kappa=0;
  bool pbc=false;
};

template<int threads>
__device__ double reduce(double value, bool maximum, double* scratch) {
  const int thread=threadIdx.x;
  if constexpr(threads==neighborThreads) {
    // a single warp of 32 threads (half a wavefront on AMD GPUs)
    for(int offset=threads/2; offset; offset/=2) {
      const double other=__shfl_down_sync(warpMask(),value,offset,32);
      value=maximum ? fmax(value,other) : value+other;
    }
    return __shfl_sync(warpMask(),value,0,32);
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

__device__ double distance(const Data& data, int center, int assigned,
                           double* displacement) {
  double squared=0;
  for(int axis=0; axis<3; ++axis) {
    displacement[axis]=data.positions[3*(data.numCenters+assigned)+axis]-
                       data.positions[3*center+axis];
    if(data.pbc) {
      // Match Pbc::distance, including its choice of image at a half-box tie.
      const double scaled=__dmul_rn(displacement[axis],data.inverseBox[axis]);
      displacement[axis]=Tools::pbc(scaled)*data.box[axis];
    }
    squared+=displacement[axis]*displacement[axis];
  }
  return sqrt(squared);
}

// The maximum score is rounded to FP64 before the reduction. Round each
// recomputed score the same way so its shift is exactly zero at the maximum,
// even when a large finite KAPPA would amplify an FMA rounding residual.
__device__ double assignmentScore(double kappa, double length) {
  return __dmul_rn(-kappa,length);
}

enum PairFailure : unsigned { nonfiniteDistance, zeroDistance, nonfiniteScore, failureKinds };
constexpr unsigned long long noFailure=ULLONG_MAX;
constexpr unsigned long long emptyList=noFailure-1;

// Lower keys win, so every launch reports the same offending pair. The two
// signed-int selection sizes leave room for the reason code and sentinels.
__device__ void recordFailure(const Data& data, int center, int assigned, PairFailure reason) {
  const auto pair=static_cast<unsigned long long>(center)*data.numAssigned+assigned;
  atomicMin(data.invalid,pair*failureKinds+reason);
}

template<int threads, bool useNeighbors>
__global__ void forwardKernel(Data data) {
  __shared__ double scratch[threads];
  const int assigned=blockIdx.x;
  const int thread=threadIdx.x;
  const int begin=useNeighbors ? data.offsets[assigned] : 0;
  const int end=useNeighbors ? data.offsets[assigned+1] : data.numCenters;
  double maximum=-INFINITY;
  for(int entry=begin+thread; entry<end; entry+=threads) {
    const int center=useNeighbors ? data.candidateCenters[entry] : entry;
    double displacement[3];
    const double length=distance(data,center,assigned,displacement);
    if(!isfinite(length) || length<=2.220446049250313e-16) {
      recordFailure(data,center,assigned,isfinite(length) ? zeroDistance : nonfiniteDistance);
    }
    if(!useNeighbors || length<=data.cutoff) {
      const double score=assignmentScore(data.kappa,length);
      if(!isfinite(score)) {
        recordFailure(data,center,assigned,nonfiniteScore);
      }
      maximum=fmax(maximum,score);
    }
  }
  maximum=reduce<threads>(maximum,true,scratch);
  if(useNeighbors && maximum==-INFINITY) {
    if(thread==0) {
      atomicMin(data.invalid,emptyList);
    }
    return;
  }
  double normalization=0;
  for(int entry=begin+thread; entry<end; entry+=threads) {
    const int center=useNeighbors ? data.candidateCenters[entry] : entry;
    double displacement[3];
    const double length=distance(data,center,assigned,displacement);
    if(!useNeighbors || length<=data.cutoff) {
      normalization+=exp(assignmentScore(data.kappa,length)-maximum);
    }
  }
  normalization=reduce<threads>(normalization,false,scratch);
  if(thread==0) {
    data.maximum[assigned]=maximum;
    data.normalization[assigned]=normalization;
  }
  for(int entry=begin+thread; entry<end; entry+=threads) {
    const int center=useNeighbors ? data.candidateCenters[entry] : entry;
    double displacement[3];
    const double length=distance(data,center,assigned,displacement);
    if(useNeighbors && length>data.cutoff) {
      continue;
    }
    const double weight=exp(assignmentScore(data.kappa,length)-maximum)/normalization;
    atomicAdd(data.occupancy+center,weight);
  }
}

template<int threads, bool useNeighbors>
__global__ void backwardKernel(Data data) {
  __shared__ double scratch[threads];
  const int assigned=blockIdx.x;
  const int thread=threadIdx.x;
  const int begin=useNeighbors ? data.offsets[assigned] : 0;
  const int end=useNeighbors ? data.offsets[assigned+1] : data.numCenters;
  const double maximum=data.maximum[assigned];
  const double normalization=data.normalization[assigned];
  double mean=0;
  for(int entry=begin+thread; entry<end; entry+=threads) {
    const int center=useNeighbors ? data.candidateCenters[entry] : entry;
    double displacement[3];
    const double length=distance(data,center,assigned,displacement);
    if(useNeighbors && length>data.cutoff) {
      continue;
    }
    const double weight=exp(assignmentScore(data.kappa,length)-maximum)/normalization;
    mean+=weight*data.sensitivity[center];
  }
  mean=reduce<threads>(mean,false,scratch);
  double assignedDerivative[3]= {0,0,0};
  double virial[9]= {};
  for(int entry=begin+thread; entry<end; entry+=threads) {
    const int center=useNeighbors ? data.candidateCenters[entry] : entry;
    double displacement[3];
    const double length=distance(data,center,assigned,displacement);
    if(useNeighbors && length>data.cutoff) {
      continue;
    }
    const double weight=exp(assignmentScore(data.kappa,length)-maximum)/normalization;
    const double radial=-data.kappa*weight*(data.sensitivity[center]-mean)/length;
    for(int axis=0; axis<3; ++axis) {
      const double gradient=radial*displacement[axis];
      atomicAdd(data.derivatives+3*center+axis,-gradient);
      assignedDerivative[axis]+=gradient;
      for(int component=0; component<3; ++component) {
        virial[3*axis+component]-=gradient*displacement[component];
      }
    }
  }
  for(int axis=0; axis<3; ++axis) {
    const double value=reduce<threads>(assignedDerivative[axis],false,scratch);
    if(thread==0) {
      data.derivatives[3*(data.numCenters+assigned)+axis]=value;
    }
  }
  for(int component=0; component<9; ++component) {
    const double value=reduce<threads>(virial[component],false,scratch);
    if(thread==0) {
      atomicAdd(data.virial+component,value);
    }
  }
}

class ReactiveVoronoiCuda final : public colvar::ReactiveVoronoiBackend {
  /// Declared first, so that it is destroyed last: the destructor enters it, so
  /// that the buffers are freed on their GPU, and the GPU of the caller is given
  /// back only after that
  CUDAHELPERS::deviceGuard destructorGuard_;
  int device_=-1;
  int numCenters_=0;
  int numAssigned_=0;
  cudaStream_t stream_=nullptr;
  deviceBuffer<double> positions_;
  deviceBuffer<double> occupancy_;
  deviceBuffer<double> maximum_;
  deviceBuffer<double> normalization_;
  deviceBuffer<double> sensitivity_;
  deviceBuffer<double> derivatives_;
  deviceBuffer<double> virial_;
  deviceBuffer<unsigned long long> invalid_;
  // pinned host copies of occupancy and invalid, written by the queued forward pass
  // (a copy to pageable memory would wait for the kernel)
  pinnedBuffer<double> hostOccupancy_;
  pinnedBuffer<unsigned long long> hostInvalid_;
  // the candidate list of NLIST: the centers of each assigned atom
  deviceBuffer<int> offsets_;
  deviceBuffer<int> candidateCenters_;
  bool haveNeighbors_=false;
  unsigned long uploadedNeighborVersion_=0;
  // the parameters of the last forward pass, used also by backward
  Data data_;
  bool useNeighbors_=false;

  void allocate() {
    CUDAHELPERS::deviceGuard guard(device_);
    cudaCheck(cudaStreamCreateWithFlags(&stream_,cudaStreamNonBlocking),"cudaStreamCreateWithFlags");
    const size_t numAtoms=static_cast<size_t>(numCenters_)+numAssigned_;
    positions_.resize(3*numAtoms);
    occupancy_.resize(numCenters_);
    maximum_.resize(numAssigned_);
    normalization_.resize(numAssigned_);
    sensitivity_.resize(numCenters_);
    derivatives_.resize(3*numAtoms);
    virial_.resize(9);
    invalid_.resize(1);
    hostOccupancy_.resize(numCenters_);
    hostInvalid_.resize(1);
  }

  // queues the forward pass, finishForward waits for it
  void queueForward(const double* positions, const Pbc& box, bool pbc, double kappa) {
    data_.numCenters=numCenters_;
    data_.numAssigned=numAssigned_;
    data_.positions=positions_.data();
    data_.occupancy=occupancy_.data();
    data_.maximum=maximum_.data();
    data_.normalization=normalization_.data();
    data_.sensitivity=sensitivity_.data();
    data_.derivatives=derivatives_.data();
    data_.virial=virial_.data();
    data_.invalid=invalid_.data();
    data_.offsets=offsets_.data();
    data_.candidateCenters=candidateCenters_.data();
    data_.pbc=pbc;
    data_.kappa=kappa;
    for(int axis=0; axis<3; ++axis) {
      data_.box[axis]=box.getBox()(axis,axis);
      data_.inverseBox[axis]=box.getInvBox()(axis,axis);
    }
    const size_t numAtoms=static_cast<size_t>(numCenters_)+numAssigned_;
    cudaCheck(cudaMemcpyAsync(positions_.data(),positions,3*numAtoms*sizeof(double),
                              cudaMemcpyHostToDevice,stream_),"cudaMemcpyAsync");
    cudaCheck(cudaMemsetAsync(occupancy_.data(),0,numCenters_*sizeof(double),stream_),"cudaMemsetAsync");
    cudaCheck(cudaMemsetAsync(invalid_.data(),0xff,sizeof(unsigned long long),stream_),"cudaMemsetAsync");
    if(useNeighbors_) {
      forwardKernel<neighborThreads,true><<<numAssigned_,neighborThreads,0,stream_>>>(data_);
    } else {
      forwardKernel<denseThreads,false><<<numAssigned_,denseThreads,0,stream_>>>(data_);
    }
    cudaCheck(cudaGetLastError(),"forwardKernel");
    cudaCheck(cudaMemcpyAsync(hostOccupancy_.data(),occupancy_.data(),numCenters_*sizeof(double),
                              cudaMemcpyDeviceToHost,stream_),"cudaMemcpyAsync");
    cudaCheck(cudaMemcpyAsync(hostInvalid_.data(),invalid_.data(),sizeof(unsigned long long),
                              cudaMemcpyDeviceToHost,stream_),"cudaMemcpyAsync");
  }

public:
  ReactiveVoronoiCuda(int centers, int assigned, int mdDevice):
    numCenters_(centers),
    numAssigned_(assigned) {
    if(centers<=0 || assigned<=0) {
      throw std::runtime_error("the GPU backend needs at least one CENTER and one ASSIGNED atom");
    }
    device_=CUDAHELPERS::chooseGPU(mdDevice).id;
    if(device_<0) {
      throw std::runtime_error("No CUDA device is available");
    }
    try {
      allocate();
    } catch(...) {
      // the members allocated so far are freed on their GPU (see destructorGuard_)
      destructorGuard_.enter(device_);
      if(stream_) {
        cudaStreamDestroy(stream_);
      }
      throw;
    }
  }

  ~ReactiveVoronoiCuda() override {
    destructorGuard_.enter(device_);
    // wait for any pending work before releasing the resources
    cudaStreamSynchronize(stream_);
    cudaStreamDestroy(stream_);
  }

  ReactiveVoronoiCuda(const ReactiveVoronoiCuda&)=delete;
  ReactiveVoronoiCuda& operator=(const ReactiveVoronoiCuda&)=delete;

  int device() const override {
    return device_;
  }

  void forward(const double* positions, const Pbc& box, bool pbc, double kappa) override {
    CUDAHELPERS::deviceGuard guard(device_);
    useNeighbors_=false;
    queueForward(positions,box,pbc,kappa);
  }

  void forwardNeighbors(const double* positions, const Pbc& box, bool pbc, double kappa,
                        double cutoff, const std::vector<int>& offsets, const std::vector<int>& centers,
                        unsigned long version) override {
    CUDAHELPERS::deviceGuard guard(device_);
    if(!haveNeighbors_ || uploadedNeighborVersion_!=version) {
      plumed_assert(offsets.size()==static_cast<size_t>(numAssigned_)+1);
      offsets_.resize(offsets.size());
      // at least one element, so that the kernels never get a null list
      candidateCenters_.resize(std::max<size_t>(centers.size(),1));
      cudaCheck(cudaMemcpyAsync(offsets_.data(),offsets.data(),offsets.size()*sizeof(int),
                                cudaMemcpyHostToDevice,stream_),"cudaMemcpyAsync");
      if(!centers.empty()) {
        cudaCheck(cudaMemcpyAsync(candidateCenters_.data(),centers.data(),centers.size()*sizeof(int),
                                  cudaMemcpyHostToDevice,stream_),"cudaMemcpyAsync");
      }
      haveNeighbors_=true;
      uploadedNeighborVersion_=version;
    }
    useNeighbors_=true;
    data_.cutoff=cutoff;
    queueForward(positions,box,pbc,kappa);
  }

  void finishForward(double* occupancy) override {
    CUDAHELPERS::deviceGuard guard(device_);
    cudaCheck(cudaStreamSynchronize(stream_),"cudaStreamSynchronize");
    std::copy(hostOccupancy_.begin(),hostOccupancy_.end(),occupancy);
    const unsigned long long invalid=hostInvalid_[0];
    if(invalid==noFailure) {
      return;
    }
    if(invalid==emptyList) {
      throw std::runtime_error("an ASSIGNED atom has no CENTER inside the candidate list");
    }
    const auto pair=invalid/failureKinds;
    const char* reasons[]= {"non-finite CENTER-ASSIGNED distance", "zero or near-zero CENTER-ASSIGNED distance",
                            "non-finite KAPPA times CENTER-ASSIGNED distance"
                           };
    throw colvar::ReactiveVoronoiBackend::InvalidPair(pair/numAssigned_,pair%numAssigned_,
        reasons[invalid%failureKinds]);
  }

  void backward(const double* sensitivity, double* derivatives, double* virial) override {
    CUDAHELPERS::deviceGuard guard(device_);
    const size_t numAtoms=static_cast<size_t>(numCenters_)+numAssigned_;
    cudaCheck(cudaMemcpyAsync(sensitivity_.data(),sensitivity,numCenters_*sizeof(double),
                              cudaMemcpyHostToDevice,stream_),"cudaMemcpyAsync");
    cudaCheck(cudaMemsetAsync(derivatives_.data(),0,3*numAtoms*sizeof(double),stream_),"cudaMemsetAsync");
    cudaCheck(cudaMemsetAsync(virial_.data(),0,9*sizeof(double),stream_),"cudaMemsetAsync");
    if(useNeighbors_) {
      backwardKernel<neighborThreads,true><<<numAssigned_,neighborThreads,0,stream_>>>(data_);
    } else {
      backwardKernel<denseThreads,false><<<numAssigned_,denseThreads,0,stream_>>>(data_);
    }
    cudaCheck(cudaGetLastError(),"backwardKernel");
    cudaCheck(cudaMemcpyAsync(derivatives,derivatives_.data(),3*numAtoms*sizeof(double),
                              cudaMemcpyDeviceToHost,stream_),"cudaMemcpyAsync");
    cudaCheck(cudaMemcpyAsync(virial,virial_.data(),9*sizeof(double),
                              cudaMemcpyDeviceToHost,stream_),"cudaMemcpyAsync");
    cudaCheck(cudaStreamSynchronize(stream_),"cudaStreamSynchronize");
  }
};

colvar::ReactiveVoronoiBackendRegistration<ReactiveVoronoiCuda> registration("CUDA");

} // namespace
} // namespace gpu_voronoi
} // namespace PLMD
