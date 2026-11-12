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
#ifndef __PLUMED_openaccPTM_ACCDevice_h
#define __PLUMED_openaccPTM_ACCDevice_h

#include <openacc.h>

#if defined(PLUMED_OPENACC_GPU) && defined(__linux__)
#include <dlfcn.h>
#endif

/*
GPU selection in the openaccPTM plugin.

The plugin runs in the process of the MD code, which may itself be using a GPU
(e.g. GROMACS with its data on the GPU). Two rules keep the two from interfering:

1. Each action (or SAXS backend) works on one GPU, chosen the first time it
   needs it and then fixed (see Device::resolve): the GPU assigned to PLUMED by
   the MD code (cmd "setGpuDeviceId"), otherwise the GPU the MD code is using
   (the one of the CUDA context current in the calling thread), otherwise the
   current OpenACC device (by default the first GPU, see ACC_DEVICE_NUM and
   CUDA_VISIBLE_DEVICES).

2. Every piece of code that uses OpenACC (kernels, transfers, allocations, and
   the destructors that free device memory) runs inside a DeviceScope for that
   GPU. The scope makes the GPU and its context current, and gives back to the
   caller its own CUDA context when it ends: what an action leaves selected
   cannot affect the MD code or the next action.
*/

namespace PLMD {
namespace ACC {

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
inline const CudaDriver& cudaDriver() {
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

/// The device of the CUDA context current in this thread, -1 if there is none
inline int currentCudaDevice() {
  const auto & cu=cudaDriver();
  void* ctx=nullptr;
  int dev=-1;
  if(cu.ok() && cu.getCurrent(&ctx)==0 && ctx && cu.getDevice(&dev)==0) {
    return dev;
  }
  return -1;
}

/// The GPU of an action, see rule 1 above
class Device {
  int id_=-1;
  bool resolved_=false;
public:
  /// The GPU to use (-1 without GPUs, when the kernels run on the host). It is
  /// chosen at the first call and then fixed. mdDevice is the GPU assigned to
  /// PLUMED by the MD code (PlumedMain::getGpuDeviceId), -1 if none.
  int resolve(const int mdDevice) {
    if(!resolved_) {
      resolved_=true;
      const int ndevices=acc_get_num_devices(acc_device_not_host);
      if(ndevices>0) {
        const int caller=currentCudaDevice();
        if(mdDevice>=0 && mdDevice<ndevices) {
          id_=mdDevice;
        } else if(caller>=0 && caller<ndevices) {
          id_=caller;
        } else {
          id_=acc_get_device_num(acc_device_not_host);
        }
      }
    }
    return id_;
  }
  /// The GPU chosen by resolve(), -1 if it was never called
  int id() const {
    return id_;
  }
};

/// Makes a GPU (and its context) current while in scope, and then gives back to
/// the caller its current CUDA context, see rule 2 above. The OpenACC runtime does
/// not make its context current again before each operation, so this is done here.
/// The OpenACC device number is not restored: the runtime is linked statically in
/// the plugin, so it is not shared with the caller, and every use of it is in a
/// DeviceScope that selects its own GPU. With device<0 nothing is done.
class DeviceScope {
  void* prevCtx_=nullptr;
  bool restore_=false;
public:
  explicit DeviceScope(const int device) {
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
  ~DeviceScope() {
    if(restore_) {
      cudaDriver().setCurrent(prevCtx_);
    }
  }
  DeviceScope(const DeviceScope&)=delete;
  DeviceScope& operator=(const DeviceScope&)=delete;
};

} // namespace ACC
} // namespace PLMD

#endif
