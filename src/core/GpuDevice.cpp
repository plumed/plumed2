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
#include "GpuDevice.h"
#include "Action.h"
#include "PlumedMain.h"
#include "tools/Log.h"

#include <cstdlib>
#if defined(__PLUMED_HAS_DLOPEN) && defined(__linux__)
#include <dlfcn.h>
#endif

namespace PLMD {

int GpuDevice::callerCudaDevice() {
#if defined(__PLUMED_HAS_DLOPEN) && defined(__linux__)
  typedef int (*getCurrent_t)(void**);
  typedef int (*getDevice_t)(int*);
  // only if the driver is already loaded; the handle is never released (the driver stays
  // loaded until the end of the process anyway)
  if(void* lib=dlopen("libcuda.so.1",RTLD_LAZY|RTLD_NOLOAD)) {
    const auto getCurrent=reinterpret_cast<getCurrent_t>(dlsym(lib,"cuCtxGetCurrent"));
    const auto getDevice=reinterpret_cast<getDevice_t>(dlsym(lib,"cuCtxGetDevice"));
    void* ctx=nullptr;
    int dev=-1;
    if(getCurrent && getDevice && getCurrent(&ctx)==0 && ctx && getDevice(&dev)==0) {
      return dev;
    }
  }
#endif
  return -1;
}

int GpuDevice::localRank(const char*& var) {
  // the variables of the launchers first: SLURM_LOCALID is also set (to 0) in the batch
  // step of a job, where a launcher like OpenMPI's mpirun may be used
  static const char* const vars[]= {"OMPI_COMM_WORLD_LOCAL_RANK","MV2_COMM_WORLD_LOCAL_RANK","MPI_LOCALRANKID",
                                    "PMI_LOCAL_RANK","PALS_LOCAL_RANKID","FLUX_TASK_LOCAL_ID","SLURM_LOCALID"
                                   };
  for(const char* v : vars) {
    if(const char* e=std::getenv(v)) {
      char* end=nullptr;
      const long r=std::strtol(e,&end,10);
      if(end!=e && *end=='\0' && r>=0) {
        var=v;
        return static_cast<int>(r);
      }
    }
  }
  return -1;
}

GpuDevice GpuDevice::choose(const int mdDevice, const Runtime& runtime) {
  GpuDevice d;
  const int n=runtime.ndevices;
  if(n<=0) {
    return d;
  }
  std::string ignored;
  if(mdDevice>=n) {
    ignored="; the GPU "+std::to_string(mdDevice)+" assigned by the MD code is not among the "+std::to_string(n)+" available";
  }
  const char* var=nullptr;
  int lr=-1;
  if(mdDevice>=0 && mdDevice<n) {
    d.id=mdDevice;
    d.why="assigned by the MD code";
  } else if(runtime.callerDevice>=0 && runtime.callerDevice<n) {
    d.id=runtime.callerDevice;
    d.why="the GPU the MD code is using";
  } else if(runtime.runtimeDeviceSetByUser && runtime.runtimeDevice>=0 && runtime.runtimeDevice<n) {
    d.id=runtime.runtimeDevice;
    d.why="chosen in the environment";
  } else if(n>1 && (lr=localRank(var))>=0) {
    d.id=lr%n;
    d.why="process "+std::to_string(lr)+" of the node, from "+var+", with "+std::to_string(n)+" GPUs";
  } else {
    d.id=(runtime.runtimeDevice>=0 && runtime.runtimeDevice<n) ? runtime.runtimeDevice : 0;
    d.why="the default";
  }
  d.why+=ignored;
  return d;
}

int GpuDevice::choose(const Action& action, const Runtime& runtime) {
  const GpuDevice d=choose(action.plumed.getGpuDeviceId(),runtime);
  if(d.id>=0) {
    action.log.printf("  %s runs on GPU %d (%s)\n",action.getLabel().c_str(),d.id,d.why.c_str());
  }
  return d.id;
}

}
