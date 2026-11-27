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
#ifndef __PLUMED_core_GpuDevice_h
#define __PLUMED_core_GpuDevice_h

#include <string>

namespace PLMD {

class Action;

/**
The GPU used by the actions that run on a GPU.

PLUMED itself does not use GPUs: the actions that do are in plugins and modules compiled with
another toolchain (plugins/openaccPTM, src/gpu_cuda_hip), each with its own runtime (OpenACC,
CUDA or HIP). They all choose the GPU here,
so that they make the same choice and write it in the log in the same way. The plugin
tells what its runtime sees, and gets the number of the GPU to use, the first that applies of:

1. the GPU assigned to PLUMED by the MD code (cmd "setGpuDeviceId"), if it is one of the
   GPUs the runtime sees;
2. the GPU the MD code is using in the thread that calls PLUMED (for CUDA, the device of the
   current context), when the plugin can tell;
3. the GPU chosen by the user for the runtime (e.g. with ACC_DEVICE_NUM);
4. with more than one GPU, and the rank of the process within its node known from the MPI
   launcher, the GPU number rank%ngpus: the processes of a node share its GPUs in turn
   (e.g. plumed driver with MPI);
5. the GPU the runtime uses by default.

GPUs are numbered as the runtime does: with CUDA_VISIBLE_DEVICES all the runtimes see the same GPUs.
*/
struct GpuDevice {
/// What the runtime of a plugin sees
  struct Runtime {
/// The number of GPUs
    int ndevices=0;
/// The GPU of the caller (rule 2), -1 if there is none or it is not known
    int callerDevice=-1;
/// The GPU the runtime would use (rules 3 and 5)
    int runtimeDevice=0;
/// The user chose runtimeDevice (rule 3)
    bool runtimeDeviceSetByUser=false;
  };
/// The number of the GPU to use, or -1 if there are no GPUs
  int id=-1;
/// How it was chosen, for the log
  std::string why;

/// The choice, given the GPU assigned by the MD code (negative if none)
  static GpuDevice choose(int mdDevice, const Runtime& runtime);
/// The choice for an action, written in its log ("<label> runs on GPU <n> (<why>)").
/// Without GPUs it returns -1 and writes nothing: the caller decides what that means.
  static int choose(const Action& action, const Runtime& runtime);
/// The device of the CUDA context that is current in this thread, i.e. the GPU the code that
/// calls PLUMED is using (Runtime::callerDevice); -1 if there is none, or without CUDA.
/// It asks the CUDA driver, if the library is already loaded in the process (PLUMED does
/// not depend on it): the plugins have their own copy of the runtimes, that would not know.
  static int callerCudaDevice();
/// The rank of this process among the processes of its node, as set in the environment
/// by the MPI launcher (-1 if it is not set); var is the variable it was read from
  static int localRank(const char*& var);
};

}

#endif
