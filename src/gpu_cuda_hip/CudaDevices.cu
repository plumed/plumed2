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
#include "cltools/CLTool.h"
#include "core/CLToolRegister.h"
#include "core/GpuDevice.h"

#include "cuda_to_hip.cuh"

#include <cstdio>

namespace PLMD {
namespace gpu_cuda_hip {

//+PLUMEDOC TOOLS cuda-devices
/*
List the GPUs that the actions of the gpu_cuda_hip module can use

This tool lists the GPUs seen by the CUDA runtime (the HIP runtime for AMD GPUs) of the
[gpu_cuda_hip module](module_gpu_cuda_hip.md), and the GPU that its actions would use without an MD code
(e.g. with [driver](driver.md)): the rules are explained in the page of the module.
It returns 0 when there is at least one GPU, 1 otherwise, so that it can be used to check whether
a GPU is available, e.g. in the regression tests (`plumed_needs_device=cuda`).

```plumed
plumed cuda-devices
```

The GPUs are only queried: the current GPU of the runtime is not changed.

*/
//+ENDPLUMEDOC

class CudaDevices :
  public CLTool {
public:
  static void registerKeywords (Keywords &keys);
  explicit CudaDevices (const CLToolOptions &co);
  int main (FILE *in, FILE *out, Communicator &pc) override;
  std::string description() const override {
    return "list the GPUs that the gpu_cuda_hip module can use";
  }
};

PLUMED_REGISTER_CLTOOL (CudaDevices, "cuda-devices")

void CudaDevices::registerKeywords (Keywords &keys) {
  CLTool::registerKeywords (keys);
  keys.addFlag ("-q", false, "print nothing: only the return value tells whether there is a GPU");
}

CudaDevices::CudaDevices (const CLToolOptions &co) :
  CLTool (co) {
  inputdata = inputType::commandline;
}

int CudaDevices::main (FILE *, FILE *out, Communicator &) {
  bool quiet;
  parseFlag ("-q", quiet);
#ifdef __HIP__
  const char *runtimeName = "HIP";
#else
  const char *runtimeName = "CUDA";
#endif
  GpuDevice::Runtime runtime;
  const cudaError_t err = cudaGetDeviceCount (&runtime.ndevices);
  if (err != cudaSuccess || runtime.ndevices <= 0) {
    if (!quiet) {
      std::fprintf (out, "no GPU available (%s runtime: %s)\n", runtimeName,
                    err != cudaSuccess ? cudaGetErrorString (err) : "no device");
    }
    return 1;
  }
  if (quiet) {
    return 0;
  }
  if (cudaGetDevice (&runtime.runtimeDevice) != cudaSuccess) {
    runtime.runtimeDevice = 0;
  }
  std::fprintf (out, "%d GPU%s seen by the %s runtime:\n", runtime.ndevices,
                runtime.ndevices > 1 ? "s" : "", runtimeName);
  for (int i = 0; i < runtime.ndevices; ++i) {
    cudaDeviceProp prop;
    if (cudaGetDeviceProperties (&prop, i) == cudaSuccess) {
      const double gib = static_cast<double> (prop.totalGlobalMem) / (1024.0 * 1024.0 * 1024.0);
#ifdef __HIP__
      std::fprintf (out, "  GPU %d: %s, %s, %.1f GiB\n", i, prop.name, prop.gcnArchName, gib);
#else
      std::fprintf (out, "  GPU %d: %s, compute capability %d.%d, %.1f GiB\n", i, prop.name,
                    prop.major, prop.minor, gib);
#endif
    } else {
      std::fprintf (out, "  GPU %d\n", i);
    }
  }
  const auto choice = GpuDevice::choose (-1, runtime);
  std::fprintf (out, "without an MD code, the actions use GPU %d (%s)\n", choice.id,
                choice.why.c_str());
  return 0;
}

} // namespace gpu_cuda_hip
} // namespace PLMD
