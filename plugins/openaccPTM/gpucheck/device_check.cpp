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
/*
Checks that the GPU plugins (openaccPTM, cudaCoord) choose their GPU with the
common rule and leave the GPU of the MD code alone (see src/ACCDevice.h here and
src/cudaHelpers.cuh in cudaCoord). It needs a node with GPUs and is run by hand
with "make check-gpu" in either plugin (not by the CI).

The program plays the MD code: it selects a GPU with the CUDA runtime and keeps a
buffer on it, as GROMACS does, then runs PLUMED with the actions of the plugins it
is given, biased so that the forces are computed too:
- openaccPTM (plumedOpenACC.*): a DISTANCE (task manager, runAllTasks and
  applyForces) and a SAXS with CHI2 (backend compute and contract), with USEGPU;
- cudaCoord (CudaCoordination.*): a CUDACOORDINATION.
After every call to PLUMED, and after PLUMED is destroyed, the GPU of the caller
must still be the current one and its buffer must still be usable.

  device_check follow plugin...   caller on the last GPU, no setGpuDeviceId:
                                  the actions must run on the caller's GPU
  device_check md plugin...       caller on the last GPU, setGpuDeviceId 0:
                                  the actions must run on GPU 0 (needs 2 GPUs)

The GPU of each action is read from the PLUMED log ("<label> runs on GPU <n>").
Exit code 0 only if every check passes.
*/
#include "plumed/wrapper/Plumed.h"

#include <cuda_runtime.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

int failures=0;

void check(const bool ok, const std::string& what) {
  std::printf("%s %s\n", ok ? "OK  " : "FAIL", what.c_str());
  if(!ok) {
    ++failures;
  }
}

/// The caller's GPU must still be current, and its buffer usable
void checkCaller(const int device, void* buffer, const std::size_t bytes, const std::string& when) {
  int current=-1;
  const cudaError_t e1=cudaGetDevice(&current);
  check(e1==cudaSuccess && current==device,
        "caller GPU " + std::to_string(device) + " is current " + when + " (found " + std::to_string(current) + ")");
  const cudaError_t e2=cudaMemset(buffer, 0, bytes);
  const cudaError_t e3=cudaDeviceSynchronize();
  check(e2==cudaSuccess && e3==cudaSuccess,
        "caller buffer usable " + when + " (" + cudaGetErrorString(e2!=cudaSuccess ? e2 : e3) + ")");
}

}

int main(int argc, char** argv) {
  const std::string mode = argc>1 ? argv[1] : "";
  if(argc<3 || (mode!="follow" && mode!="md")) {
    std::fprintf(stderr, "usage: %s follow|md plugin.so [plugin.so...]\n", argv[0]);
    return 2;
  }
  const std::vector<std::string> plugins(argv+2, argv+argc);
  int ngpus=0;
  if(cudaGetDeviceCount(&ngpus)!=cudaSuccess || ngpus<1) {
    std::fprintf(stderr, "no GPU available\n");
    return 2;
  }
  if(mode=="md" && ngpus<2) {
    std::printf("SKIP md mode needs at least 2 GPUs\n");
    return 0;
  }
  // the MD code: its GPU (the last one, so that it is not the default OpenACC device) and a buffer on it
  const int callerGpu=ngpus-1;
  const int expectedGpu = mode=="md" ? 0 : callerGpu;
  cudaSetDevice(callerGpu);
  const std::size_t bytes=1<<20;
  void* buffer=nullptr;
  check(cudaMalloc(&buffer, bytes)==cudaSuccess, "caller allocates on GPU " + std::to_string(callerGpu));

  const std::string logfile="device_check_" + mode + ".log";
  std::vector<std::string> labels;
  {
    PLMD::Plumed p;
    int natoms=4;
    double timestep=0.001;
    p.cmd("setNatoms", &natoms);
    p.cmd("setMDEngine", "device_check");
    p.cmd("setTimestep", &timestep);
    p.cmd("setLogFile", logfile.c_str());
    if(mode=="md") {
      int md=0;
      p.cmd("setGpuDeviceId", &md);
    }
    p.cmd("init");
    checkCaller(callerGpu, buffer, bytes, "after init");
    for(const auto & plugin : plugins) {
      p.cmd("readInputLine", ("LOAD FILE=" + plugin).c_str());
      if(plugin.find("plumedOpenACC")!=std::string::npos) {
        // a vector of distances: a single DISTANCE is computed on the CPU also with USEGPU
        p.cmd("readInputLine", "d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 USEGPU");
        p.cmd("readInputLine", "dsum: SUM ARG=d PERIODIC=NO");
        p.cmd("readInputLine", "RESTRAINT ARG=dsum AT=0.5 KAPPA=10");
        p.cmd("readInputLine", "s: SAXS ATOMS=1-4 NOPBC USEGPU CHI2 "
              "PARAMETERS1=1.0 PARAMETERS2=2.0 PARAMETERS3=1.5 PARAMETERS4=0.5 "
              "QVALUE1=0.03 EXPINT1=8.3 EXPERR1=0.2 QVALUE2=0.12 EXPINT2=8.8 EXPERR2=0.25 "
              "QVALUE3=0.25 EXPINT3=5.2 EXPERR3=0.15 QVALUE4=0.40 EXPINT4=2.0 EXPERR4=0.10");
        p.cmd("readInputLine", "RESTRAINT ARG=s.chi2 AT=0 KAPPA=1");
        labels.push_back("d");
        labels.push_back("s");
      } else if(plugin.find("CudaCoordination")!=std::string::npos) {
        p.cmd("readInputLine", "c: CUDACOORDINATION GROUPA=1-4 R_0=0.3 D_MAX=1.0 NOPBC");
        p.cmd("readInputLine", "RESTRAINT ARG=c AT=0 KAPPA=1");
        labels.push_back("c");
      } else {
        check(false, "known plugin " + plugin);
      }
    }
    checkCaller(callerGpu, buffer, bytes, "after the input");

    std::vector<double> pos= {0,0,0, 0.3,0,0, 0,0.4,0, 0.2,0.2,0.3};
    std::vector<double> masses(natoms,1.0), charges(natoms,0.0), forces(3*natoms);
    std::vector<double> box(9,0.0), virial(9);
    for(int step=0; step<3; ++step) {
      std::fill(forces.begin(), forces.end(), 0.0);
      std::fill(virial.begin(), virial.end(), 0.0);
      pos[0]+=0.01;
      p.cmd("setStep", &step);
      p.cmd("setPositions", pos.data());
      p.cmd("setMasses", masses.data());
      p.cmd("setCharges", charges.data());
      p.cmd("setBox", box.data());
      p.cmd("setForces", forces.data());
      p.cmd("setVirial", virial.data());
      p.cmd("calc");
      checkCaller(callerGpu, buffer, bytes, "after step " + std::to_string(step));
      bool finite=true, nonzero=false;
      for(const auto f : forces) {
        finite = finite && f==f;
        nonzero = nonzero || f!=0.0;
      }
      check(finite && nonzero, "forces finite and non zero at step " + std::to_string(step));
    }
  }
  checkCaller(callerGpu, buffer, bytes, "after PLUMED is destroyed");

  std::ifstream log(logfile);
  std::stringstream text;
  text << log.rdbuf();
  for(const auto & label : labels) {
    const std::string expected=label + " runs on GPU " + std::to_string(expectedGpu);
    check(text.str().find(expected)!=std::string::npos, "log says \"" + expected + "\"");
  }

  cudaFree(buffer);
  std::printf("%s: %d failure(s)\n", mode.c_str(), failures);
  return failures==0 ? 0 : 1;
}
