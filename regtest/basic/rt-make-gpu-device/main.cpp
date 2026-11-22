#include "plumed/core/GpuDevice.h"
#include "plumed/core/ActionAnyorder.h"
#include "plumed/core/ActionRegister.h"
#include "plumed/core/PlumedMain.h"
#include <cstdlib>
#include <fstream>
#include <string>
#if defined(__PLUMED_HAS_DLOPEN) && defined(__linux__)
#include <dlfcn.h>
#endif

using namespace PLMD;

// An action that chooses its GPU as the plugins do, with an invented runtime
class GpuDeviceTest : public ActionAnyorder {
public:
  static void registerKeywords(Keywords& keys) {
    ActionAnyorder::registerKeywords(keys);
    keys.add("compulsory","NDEVICES","the number of GPUs the runtime sees");
    keys.add("compulsory","CALLER","-1","the GPU of the caller");
  }
  explicit GpuDeviceTest(const ActionOptions& ao): Action(ao), ActionAnyorder(ao) {
    GpuDevice::Runtime r;
    parse("NDEVICES",r.ndevices);
    parse("CALLER",r.callerDevice);
    checkRead();
    log.printf("  chosen: %d\n",GpuDevice::choose(*this,r));
  }
};
PLUMED_REGISTER_ACTION(GpuDeviceTest,"GPU_DEVICE_TEST")

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
#if defined(__PLUMED_HAS_DLOPEN) && defined(__linux__)
  // a fake CUDA driver (libcuda.so.1, built by the config of this test) with a current
  // context on device 3, and then none; without dlopen the reference lines are written here
  if(void* lib=dlopen("./libcuda.so.1",RTLD_NOW|RTLD_GLOBAL)) {
    auto setDevice=reinterpret_cast<void(*)(int)>(dlsym(lib,"fakeSetDevice"));
    setDevice(3);
    out<<"caller CUDA device, fake driver: "<<GpuDevice::callerCudaDevice()<<"\n";
    setDevice(-1);
    out<<"caller CUDA device, fake driver, no context: "<<GpuDevice::callerCudaDevice()<<"\n";
  } else {
    out<<"no fake driver: "<<dlerror()<<"\n";
  }
#else
  out<<"caller CUDA device, fake driver: 3\n";
  out<<"caller CUDA device, fake driver, no context: -1\n";
#endif
  // the choice of an action, with the GPU assigned by the MD code, written in its log
  {
    PlumedMain pm;
    int natoms=1;
    int md=2;
    pm.cmd("setNatoms",&natoms);
    pm.cmd("setLogFile","gpu.log");
    pm.cmd("setGpuDeviceId",&md);
    pm.cmd("init");
    pm.cmd("readInputLine","a: GPU_DEVICE_TEST NDEVICES=4");
    pm.cmd("readInputLine","b: GPU_DEVICE_TEST NDEVICES=2 CALLER=1");
    pm.cmd("readInputLine","c: GPU_DEVICE_TEST NDEVICES=0");
  }
  std::ifstream log("gpu.log");
  for(std::string line; std::getline(log,line);) {
    if(line.find("runs on GPU")!=std::string::npos || line.find("chosen:")!=std::string::npos) {
      out<<line<<"\n";
    }
  }
  return 0;
}
