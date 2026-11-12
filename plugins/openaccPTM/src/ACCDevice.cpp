/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ACCDevice.h"

#include "plumed/core/Action.h"
#include "plumed/core/PlumedMain.h"
#include "plumed/tools/Log.h"

#include <openacc.h>

#if defined(PLUMED_OPENACC_GPU) && defined(__linux__)
#include <dlfcn.h>
#endif

namespace PLMD {
namespace ACC {

namespace {

/// The few functions of the CUDA driver needed to read and restore the CUDA
/// context of the caller. They are taken from the driver library that is already
/// loaded (by the caller or by the OpenACC runtime), so that the plugin does not
/// depend on it at link time. Only with GPUs (PLUMED_OPENACC_GPU, set by the Makefile).
struct CudaDriver {
  typedef int (*getCurrent_t)(void**);
  typedef int (*setCurrent_t)(void*);
  typedef int (*getDevice_t)(int*);
  getCurrent_t getCurrent=nullptr;
  setCurrent_t setCurrent=nullptr;
  getDevice_t getDevice=nullptr;
  bool ok() const {
    return getCurrent && setCurrent && getDevice;
  }
};

/// The CUDA driver, or an empty one if it is not (yet) loaded or not available
const CudaDriver& cudaDriver() {
  static CudaDriver driver;
#if defined(PLUMED_OPENACC_GPU) && defined(__linux__)
  if(!driver.ok()) {
    // never released: the driver stays loaded until the end of the process anyway
    if(void* lib=dlopen("libcuda.so.1",RTLD_LAZY|RTLD_NOLOAD)) {
      CudaDriver d;
      d.getCurrent=reinterpret_cast<CudaDriver::getCurrent_t>(dlsym(lib,"cuCtxGetCurrent"));
      d.setCurrent=reinterpret_cast<CudaDriver::setCurrent_t>(dlsym(lib,"cuCtxSetCurrent"));
      d.getDevice=reinterpret_cast<CudaDriver::getDevice_t>(dlsym(lib,"cuCtxGetDevice"));
      if(d.ok()) {
        driver=d;
      }
    }
  }
#endif
  return driver;
}

}

int currentCudaDevice() {
  const auto & cu=cudaDriver();
  void* ctx=nullptr;
  int dev=-1;
  if(cu.ok() && cu.getCurrent(&ctx)==0 && ctx && cu.getDevice(&dev)==0) {
    return dev;
  }
  return -1;
}

int Device::resolve(const Action& action) {
  if(!resolved_) {
    resolved_=true;
    const int ndevices=acc_get_num_devices(acc_device_not_host);
    if(ndevices>0) {
      const int md=action.plumed.getGpuDeviceId();
      const int caller=currentCudaDevice();
      if(md>=0 && md<ndevices) {
        id_=md;
      } else if(caller>=0 && caller<ndevices) {
        id_=caller;
      } else {
        id_=acc_get_device_num(acc_device_not_host);
      }
      action.log.printf("  %s runs on GPU %d\n",action.getLabel().c_str(),id_);
    } else {
      action.log.printf("  %s runs on the host (no GPU available to OpenACC)\n",action.getLabel().c_str());
    }
  }
  return id_;
}

DeviceScope::DeviceScope(const int device) {
  if(device<0) {
    return;
  }
  const auto & cu=cudaDriver();
  restore_ = cu.ok() && cu.getCurrent(&prevCtx_)==0;
  if(acc_get_device_num(acc_device_not_host)!=device) {
    acc_set_device_num(device,acc_device_not_host);
  }
  if(restore_) {
    // null if the runtime has not initialized the device yet: it makes its
    // context current when it does
    void* ctx=acc_get_current_cuda_context();
    void* cur=nullptr;
    if(ctx && cu.getCurrent(&cur)==0 && cur!=ctx) {
      cu.setCurrent(ctx);
    }
  }
}

DeviceScope::~DeviceScope() {
  if(restore_) {
    cudaDriver().setCurrent(prevCtx_);
  }
}

} // namespace ACC
} // namespace PLMD
