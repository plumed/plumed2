#include "plumed/core/GpuDevice.h"
#include <cstdlib>
#include <fstream>

using namespace PLMD;

// The choice of the GPU made for the plugins (GpuDevice::choose), without GPUs:
// the plugin tells what its runtime sees, here it is invented.

static const char* const rankVars[]= {"OMPI_COMM_WORLD_LOCAL_RANK","MV2_COMM_WORLD_LOCAL_RANK","MPI_LOCALRANKID",
                                      "PMI_LOCAL_RANK","PALS_LOCAL_RANKID","FLUX_TASK_LOCAL_ID","SLURM_LOCALID"
                                     };

static void clearRank() {
  for(const char* v : rankVars) {
    unsetenv(v);
  }
}

int main() {
  std::ofstream out("output");
  auto test=[&](const char* what, int md, int ndevices, int caller, int runtimeDevice, bool byUser) {
    GpuDevice::Runtime r;
    r.ndevices=ndevices;
    r.callerDevice=caller;
    r.runtimeDevice=runtimeDevice;
    r.runtimeDeviceSetByUser=byUser;
    const GpuDevice d=GpuDevice::choose(md,r);
    out<<what<<": GPU "<<d.id<<" ("<<d.why<<")\n";
  };
  clearRank();
  test("no GPUs",2,0,1,0,false);
  test("MD code",2,4,1,0,true);
  test("MD code, not available",5,4,1,0,false);
  test("MD code, not available, no caller",5,4,-1,3,false);
  test("caller",-1,4,1,0,true);
  test("caller, not available",-1,4,7,2,true);
  test("environment",-1,4,-1,2,true);
  test("default",-1,4,-1,0,false);
  test("default, one GPU",-1,1,-1,0,false);
  test("default, wrong runtime device",-1,2,-1,9,false);
  // the rank of the process in its node
  setenv("SLURM_LOCALID","0",1);
  setenv("OMPI_COMM_WORLD_LOCAL_RANK","5",1);
  test("rank, launcher before SLURM",-1,4,-1,0,false);
  test("rank, one GPU",-1,1,-1,0,false);
  test("rank, MD code first",3,4,-1,0,false);
  test("rank, caller first",-1,4,2,0,false);
  test("rank, environment first",-1,4,-1,3,true);
  clearRank();
  setenv("SLURM_LOCALID","3",1);
  test("rank, SLURM",-1,2,-1,0,false);
  setenv("SLURM_LOCALID","x3",1);
  test("rank, not a number",-1,2,-1,1,false);
  setenv("SLURM_LOCALID","-1",1);
  test("rank, negative",-1,2,-1,1,false);
  clearRank();
  // only that it can be called: there is no CUDA context here
  out<<"caller CUDA device: "<<GpuDevice::callerCudaDevice()<<"\n";
  return 0;
}
