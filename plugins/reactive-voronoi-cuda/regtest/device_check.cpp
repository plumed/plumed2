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
#include <plumed/wrapper/Plumed.h>
#include <cuda_runtime.h>
#include <dlfcn.h>
#include <array>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

void check(cudaError_t status) {
  if(status!=cudaSuccess) {
    throw std::runtime_error(cudaGetErrorString(status));
  }
}

void currentIs(int expected) {
  int current=-1;
  check(cudaGetDevice(&current));
  if(current!=expected) {
    throw std::runtime_error("The MD current device changed");
  }
}

std::array<double,12> evaluate(const std::string& backend, int mdDevice,
                               const std::string& log) {
  PLMD::Plumed plumed;
  int atoms=4;
  int precision=8;
  double timestep=0.0005;
  plumed.cmd("setNatoms",&atoms);
  plumed.cmd("setRealPrecision",&precision);
  plumed.cmd("setTimestep",&timestep);
  plumed.cmd("setLogFile",log.c_str());
  if(mdDevice>=0) {
    plumed.cmd("setGpuDeviceId",&mdDevice);
  }
  plumed.cmd("init");
  const std::string line="q: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3,4"
                         " REFERENCE=1 KAPPA=5 POWER=2"+backend;
  plumed.cmd("readInputLine",line.c_str());
  plumed.cmd("readInputLine","r: RESTRAINT ARG=q AT=0.3 KAPPA=3");
  std::array<double,12> positions{{0.1,0.2,0.1, 1.5,0.3,0.2, 0.4,0.3,0.2, 1.1,0.5,0.2}};
  std::array<double,12> forces{};
  std::array<double,4> masses{{1,1,1,1}};
  std::array<double,9> box{{3,0,0,0,3,0,0,0,3}},virial{};
  int step=0;
  plumed.cmd("setStep",&step);
  plumed.cmd("setPositions",positions.data());
  plumed.cmd("setMasses",masses.data());
  plumed.cmd("setBox",box.data());
  plumed.cmd("setForces",forces.data());
  plumed.cmd("setVirial",virial.data());
  plumed.cmd("calc");
  return forces;
}

int main(int argc,char** argv) {
  if(argc!=2) {
    return 2;
  }
  int count=0;
  check(cudaGetDeviceCount(&count));
  if(count<2) {
    std::cout<<"SKIP two GPUs are needed for distinct MD/CV contexts\n";
    return 77;
  }
  check(cudaSetDevice(1));
  double* mdMemory=nullptr;
  check(cudaMalloc(&mdMemory,16*sizeof(double)));
  const std::string library=argv[1];
  const auto reference=evaluate("",-1,"cpu.log");
  double maximum=0;
  for(int requested: {
        0,-1
      }) {
    const std::string log=requested==0 ? "md-selected.log" : "md-current.log";
    const auto forces=evaluate(" BACKEND=CUDA CUDA_LIBRARY="+library,requested,log);
    currentIs(1);
    check(cudaMemset(mdMemory,0,16*sizeof(double)));
    check(cudaDeviceSynchronize());
    for(unsigned index=0; index<forces.size(); ++index) {
      maximum=std::max(maximum,std::abs(forces[index]-reference[index]));
    }
    std::ifstream stream(log);
    const std::string text((std::istreambuf_iterator<char>(stream)), {});
    const std::string device=requested==0 ? "device 0" : "device 1";
    if(text.find("selected CUDA ("+device)==std::string::npos) {
      throw std::runtime_error("The logged CV device differs from the requested device");
    }
  }
  for(int repeat=0; repeat<8; ++repeat) {
    void* libraryHandle=dlopen(library.c_str(),RTLD_NOW|RTLD_LOCAL);
    if(!libraryHandle) {
      throw std::runtime_error(dlerror());
    }
    const auto create=reinterpret_cast<void*(*)(int,int,int)>(
                        dlsym(libraryHandle,"reactive_voronoi_create_v2"));
    const auto destroy=reinterpret_cast<void(*)(void*)>(
                         dlsym(libraryHandle,"reactive_voronoi_destroy_v1"));
    const auto device=reinterpret_cast<int(*)(void*)>(
                        dlsym(libraryHandle,"reactive_voronoi_device_v1"));
    const auto forward=reinterpret_cast<int(*)(void*,const double*,const double*,int,double,double*)>(
                         dlsym(libraryHandle,"reactive_voronoi_forward_v1"));
    const auto backward=reinterpret_cast<int(*)(void*,const double*,double*,double*)>(
                          dlsym(libraryHandle,"reactive_voronoi_backward_v1"));
    if(!create || !destroy || !device || !forward || !backward) {
      throw std::runtime_error("CUDA library ABI is incomplete");
    }
    if(create(2,2,count)!=nullptr) {
      throw std::runtime_error("An invalid requested GPU was accepted");
    }
    currentIs(1);
    void* context=create(2,2,0);
    if(!context || device(context)!=0) {
      throw std::runtime_error("CUDA context has the wrong owner");
    }
    currentIs(1);
    const double positions[12]= {0.1,0.2,0.1, 1.5,0.3,0.2, 0.4,0.3,0.2, 1.1,0.5,0.2};
    const double box[3]= {3,3,3};
    const double sensitivity[2]= {1.2,-0.7};
    double occupancy[2],derivatives[12],virial[9];
    if(forward(context,positions,box,1,5,occupancy)!=0 ||
        backward(context,sensitivity,derivatives,virial)!=0) {
      throw std::runtime_error("CUDA context evaluation failed");
    }
    currentIs(1);
    destroy(context);
    currentIs(1);
    dlclose(libraryHandle);
    check(cudaMemset(mdMemory,0,16*sizeof(double)));
    check(cudaDeviceSynchronize());
    currentIs(1);
  }
  check(cudaFree(mdMemory));
  if(maximum>2e-10) {
    throw std::runtime_error("Device-selection force parity failed");
  }
  std::cout<<"PASS requested/current device selection, restore and MD allocation lifetime\n"
           <<"maximum_force_difference "<<maximum<<"\n";
}
