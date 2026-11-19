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
#ifndef __PLUMED_colvar_ReactiveVoronoiCuda_h
#define __PLUMED_colvar_ReactiveVoronoiCuda_h

#include "tools/DLLoader.h"
#ifdef __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif
#include <stdexcept>
#include <string>
#include <vector>

namespace PLMD {
namespace colvar {

/// Optional runtime backend. CUDA headers are confined to the plugin.
class ReactiveVoronoiCuda {
  DLLoader loader_;
  void* library_=nullptr;
  void* state_=nullptr;
  int device_=-1;
  unsigned long uploadedNeighborVersion_=0;
  using Create=void*(*)(int,int);
  using CreateOnDevice=void*(*)(int,int,int);
  using Destroy=void(*)(void*);
  using Forward=int(*)(void*,const double*,const double*,int,double,double*);
  using ForwardNeighbors=int(*)(void*,const double*,const double*,int,double,double*,
                                double,const int*,const int*,int);
  using Backward=int(*)(void*,const double*,double*,double*);
  using Device=int(*)(void*);
  using DeviceCount=int(*)();
  Destroy destroy_=nullptr;
  Forward forward_=nullptr;
  ForwardNeighbors forwardNeighbors_=nullptr;
  Backward backward_=nullptr;

  template<class Function>
  Function symbol(const char* name, bool required=true) {
#ifdef __PLUMED_HAS_DLOPEN
    void* value=dlsym(library_,name);
    if(!value && required) {
      throw std::runtime_error(std::string("Missing CUDA backend symbol ")+name);
    }
    return reinterpret_cast<Function>(value);
#else
    (void)name;
    (void)required;
    throw std::runtime_error("Dynamic loading is unavailable in this PLUMED build");
#endif
  }

public:
  ReactiveVoronoiCuda(const std::string& path, int centers, int assigned,
                      int mdDevice=-1, bool needsNeighbors=false) {
    if(!DLLoader::installed()) {
      throw std::runtime_error("Dynamic loading is unavailable in this PLUMED build");
    }
    library_=loader_.load(path);
    const auto create=symbol<Create>("reactive_voronoi_create_v1");
    const auto createOnDevice=symbol<CreateOnDevice>("reactive_voronoi_create_v2",false);
    const auto deviceCount=symbol<DeviceCount>("reactive_voronoi_device_count_v1",false);
    const auto device=symbol<Device>("reactive_voronoi_device_v1",false);
    destroy_=symbol<Destroy>("reactive_voronoi_destroy_v1");
    forward_=symbol<Forward>("reactive_voronoi_forward_v1");
    if(needsNeighbors) {
      forwardNeighbors_=symbol<ForwardNeighbors>("reactive_voronoi_forward_nlist_v1",false);
      if(!forwardNeighbors_) {
        throw std::runtime_error("The CUDA library does not support NLIST");
      }
    }
    backward_=symbol<Backward>("reactive_voronoi_backward_v1");
    if(deviceCount && deviceCount()==0) {
      throw std::runtime_error("No CUDA device is available");
    }
    if(mdDevice>=0 && !createOnDevice) {
      throw std::runtime_error("This CUDA library needs the v2 factory to honor setGpuDeviceId");
    }
    state_=createOnDevice ? createOnDevice(centers,assigned,mdDevice) :
           create(centers,assigned);
    if(!state_) {
      throw std::runtime_error("CUDA backend context creation failed on the requested or current device");
    }
    if(device) {
      device_=device(state_);
    }
  }

  ~ReactiveVoronoiCuda() {
    if(state_) {
      destroy_(state_);
    }
  }
  ReactiveVoronoiCuda(const ReactiveVoronoiCuda&)=delete;
  ReactiveVoronoiCuda& operator=(const ReactiveVoronoiCuda&)=delete;

  int device() const {
    return device_;
  }

  void forward(const double* positions, const double* box, bool pbc,
               double kappa, double* occupancy) {
    if(forward_(state_,positions,box,pbc,kappa,occupancy)!=0) {
      throw std::runtime_error("CUDA reactive Voronoi forward failed; check geometry and CUDA diagnostics");
    }
  }

  void forwardNeighbors(const double* positions, const double* box, bool pbc,
                        double kappa, double* occupancy, double cutoff,
                        const std::vector<int>& offsets, const std::vector<int>& centers,
                        unsigned long version) {
    const bool upload=uploadedNeighborVersion_!=version;
    const int status=forwardNeighbors_(state_,positions,box,pbc,kappa,occupancy,cutoff,
                                       upload ? offsets.data() : nullptr,
                                       upload ? centers.data() : nullptr,
                                       upload ? static_cast<int>(centers.size()) : 0);
    if(status==2) {
      throw std::runtime_error("an ASSIGNED atom has no CENTER inside the candidate list");
    }
    if(status!=0) {
      throw std::runtime_error("CUDA reactive Voronoi NLIST forward failed; check geometry and CUDA diagnostics");
    }
    uploadedNeighborVersion_=version;
  }

  void backward(const double* sensitivity, double* derivatives, double* virial) {
    if(backward_(state_,sensitivity,derivatives,virial)!=0) {
      throw std::runtime_error("CUDA reactive Voronoi backward failed");
    }
  }
};
}
}
#endif
