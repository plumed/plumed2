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

/*
GPU selection in the openaccPTM plugin.

The plugin runs in the process of the MD code, which may itself be using a GPU
(e.g. GROMACS with its data on the GPU). Two rules keep the two from interfering,
and every part of the plugin follows them in the same way:

1. Each action works on one GPU, held in an ACC::Device and chosen the first
   time it is needed with Device::resolve(action), which also writes it in the
   log of the action ("<label> runs on GPU <n>"): the GPU assigned to PLUMED by
   the MD code (cmd "setGpuDeviceId"), otherwise the GPU the MD code is using
   (the one of the CUDA context current in the calling thread), otherwise the
   current OpenACC device (by default the first GPU, see ACC_DEVICE_NUM and
   CUDA_VISIBLE_DEVICES). It is then fixed.

2. Every piece of code that uses OpenACC (kernels, transfers, allocations, and
   the destructors that free device memory) runs inside a DeviceScope for that
   GPU. The scope makes the GPU and its context current, and gives back to the
   caller its own CUDA context when it ends: what an action leaves selected
   cannot affect the MD code or the next action.

The actions that use the task manager (ACCParallelTaskManager.h) get both for
free. Code that uses OpenACC outside of it (e.g. the SAXS/SANS backend in
OpenACCDebye.cpp) keeps an ACC::Device, calls resolve(action) before its first
use of the GPU, and puts every use, and the release of its device memory in
the destructor, in a DeviceScope(device.id()).
*/

namespace PLMD {

class Action;

namespace ACC {

/// The device of the CUDA context current in this thread, -1 if there is none
/// (or without GPUs)
int currentCudaDevice();

/// The GPU of an action, see rule 1 above
class Device {
  int id_=-1;
  bool resolved_=false;
public:
  /// The GPU to use for an action (-1 without GPUs, when the kernels run on the
  /// host). It is chosen at the first call, written in the log of the action,
  /// and then fixed.
  int resolve(const Action& action);
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
  explicit DeviceScope(int device);
  ~DeviceScope();
  DeviceScope(const DeviceScope&)=delete;
  DeviceScope& operator=(const DeviceScope&)=delete;
};

} // namespace ACC
} // namespace PLMD

#endif
