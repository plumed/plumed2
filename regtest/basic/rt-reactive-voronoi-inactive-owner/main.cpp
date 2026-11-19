#include "plumed/wrapper/Plumed.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void test(bool cells, bool strideOwner) {
  const int extra=cells ? 128 : 0;
  const int centers=2+extra, assigned=1+2*extra;
  const int atoms=centers+assigned;
  std::vector<double> x(3*atoms,0.0), forces(3*atoms), masses(atoms,1.0);
  const double length=cells ? 12.0 : 3.0;
  // One assigned atom sees both centers across the periodic seam.
  x[0]=0.1;
  x[3]=length-0.13;
  x[3*centers]=length-0.04;
  for(int i=0; i<extra; ++i) {
    const std::array<double,3> center{{1.0+i%8,1.0+(i/8)%8,1.0+i/64}};
    for(int k=0; k<3; ++k) {
      x[3*(i+2)+k]=center[k];
      x[3*(centers+2*i+1)+k]=center[k]+(k==0 ? 0.12 : 0.0);
      x[3*(centers+2*i+2)+k]=center[k]+(k==1 ? 0.14 : 0.0);
    }
  }
  std::array<double,9> box{{length,0,0,0,length,0,0,0,length}}, virial{};
  const std::string label=std::string(cells ? "cells" : "scan")+(strideOwner ? "-stride" : "-unused");
  PLMD::Plumed p;
  p.cmd("setRealPrecision",8);
  p.cmd("setNatoms",atoms);
  p.cmd("setTimestep",0.001);
  const std::string logfile=label+".log";
  p.cmd("setLogFile",logfile.c_str());
  p.cmd("init");
  const std::string parameters=": VORONOI_COORDINATION CENTERS=1-"+std::to_string(centers)+
                               " ASSIGNED="+std::to_string(centers+1)+"-"+std::to_string(atoms)+
                               " REFERENCE=0 POWER=2 NLIST NL_CUTOFF=0.30 NL_SKIN=0.10 NL_STRIDE=5";
  const auto read=[&](const std::string& line) {
    p.cmd("readInputLine",line.c_str());
  };
  read("owner"+parameters+" KAPPA=4");
  read("active"+parameters+" KAPPA=5");
  read("reference"+parameters+" KAPPA=5 NO_SHARE");
  if(strideOwner) {
    read("PRINT ARG=owner STRIDE=100 FILE="+label+".colvar");
  }
  read("bias: RESTRAINT ARG=active,reference AT=0.25,0.25 KAPPA=2,2");
  double active=0.0, reference=0.0;
  long rank=-1;
  p.cmd("getDataRank active",&rank);
  p.cmd("setMemoryForData active",&active);
  p.cmd("getDataRank reference",&rank);
  p.cmd("setMemoryForData reference",&reference);
  for(int step=0; step<3; ++step) {
    // Step 1 changes the box and makes the first image outside the owner's
    // stored candidate cutoff. The active action still sees both centers.
    box[0]=length-(step ? 0.42 : 0.0);
    x[0]=step ? 0.32 : 0.1;
    x[3]=box[0]-0.22;
    x[3*centers]=box[0]+0.12;
    if(step==0) {
      x[0]=0.1;
      x[3]=length-0.13;
      x[3*centers]=length-0.04;
    }
    p.cmd("setStep",step);
    p.cmd("setPositions",x.data());
    p.cmd("setMasses",masses.data());
    p.cmd("setForces",forces.data());
    p.cmd("setBox",box.data());
    p.cmd("setVirial",virial.data());
    p.cmd("calc");
    if(!std::isfinite(active) || !std::isfinite(reference) ||
        std::abs(active-reference)>1e-10) {
      throw std::runtime_error(label+" active and independent CV differ at step "+std::to_string(step));
    }
  }
}

}

int main() {
  try {
    for(bool cells: {
          false,true
        }) {
      for(bool strideOwner: {
            false,true
          }) {
        test(cells,strideOwner);
      }
    }
    std::ofstream("output")<<"PASS reactive Voronoi uses the active PBC for shared list rebuilds\n";
  } catch(const std::exception& error) {
    std::cerr<<error.what()<<'\n';
    return 1;
  }
}
