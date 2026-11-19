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

std::array<double,12> evaluate(const std::string& library, int mdDevice,
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
  if(!library.empty()) {
    plumed.cmd("readInputLine","LOAD FILE="+library);
  }
  const std::string line="q: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3,4"
                         " REFERENCE=1 KAPPA=5 POWER=2"+
                         std::string(library.empty() ? "" : " BACKEND=CUDA");
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
  if(count==0) {
    std::cout<<"SKIP no visible CUDA device\n";
    return 77;
  }
  // With two devices, verify that offloading to device 0 restores device 1.
  const int current=count>1 ? 1 : 0;
  check(cudaSetDevice(current));
  double* mdMemory=nullptr;
  check(cudaMalloc(&mdMemory,16*sizeof(double)));
  const std::string library=argv[1];
  const auto reference=evaluate("",-1,"cpu.log");
  double maximum=0;
  for(int requested: {
        0,-1
      }) {
    const std::string log=requested==0 ? "md-selected.log" : "md-current.log";
    const auto forces=evaluate(library,requested,log);
    currentIs(current);
    check(cudaMemset(mdMemory,0,16*sizeof(double)));
    check(cudaDeviceSynchronize());
    for(unsigned index=0; index<forces.size(); ++index) {
      maximum=std::max(maximum,std::abs(forces[index]-reference[index]));
    }
    std::ifstream stream(log);
    const std::string text((std::istreambuf_iterator<char>(stream)), {});
    const std::string device="device "+std::to_string(requested==0 ? 0 : current);
    if(text.find("selected CUDA ("+device)==std::string::npos) {
      throw std::runtime_error("The logged CV device differs from the requested device");
    }
  }
  try {
    evaluate(library,count,"invalid-device.log");
    throw std::runtime_error("An invalid requested GPU was accepted");
  } catch(const std::exception& error) {
    if(std::string(error.what()).find("CUDA initialization failed")==std::string::npos) {
      throw;
    }
  }
  currentIs(current);
  for(int repeat=0; repeat<8; ++repeat) {
    const auto forces=evaluate(library,0,"repeat.log");
    currentIs(current);
    for(unsigned index=0; index<forces.size(); ++index) {
      maximum=std::max(maximum,std::abs(forces[index]-reference[index]));
    }
    check(cudaMemset(mdMemory,0,16*sizeof(double)));
    check(cudaDeviceSynchronize());
  }
  check(cudaFree(mdMemory));
  if(maximum>2e-10) {
    throw std::runtime_error("Device-selection force parity failed");
  }
  std::cout<<"PASS requested/current device selection, restore and MD allocation lifetime\n"
           <<"visible_devices "<<count<<"\n"
           <<"maximum_force_difference "<<maximum<<"\n";
}
