#include "plumed/wrapper/Plumed.h"
#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

int main() {
  std::ofstream output("output");
  for(int test=0; test<4; ++test) {
    const bool energy=test!=1;
    const bool hrex=test!=0;
    PLMD::Plumed p;
    p.cmd("GREX init");
    int natoms=2;
    double timestep=0.001;
    p.cmd("setNatoms",&natoms);
    p.cmd("setTimestep",&timestep);
    const std::string logfile="log-"+std::to_string(test);
    p.cmd("setLogFile",logfile.c_str());
    p.cmd("init");
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    if(energy) {
      p.cmd("readInputLine","e: ENERGY");
      if(test!=3) {
        p.cmd("readInputLine","b: RESTRAINT ARG=e AT=0 KAPPA=2");
      }
    } else {
      p.cmd("readInputLine","b: RESTRAINT ARG=d AT=0 KAPPA=2");
    }
    std::array<double,6> positions{{0,0,0,1,0,0}},forces{};
    std::array<double,2> masses{{1,1}};
    std::array<double,9> box{{10,0,0,0,10,0,0,0,10}},virial{};
    int step=0;
    double own=2,foreign=4,bias=0;
    p.cmd("setStep",&step);
    p.cmd("setPositions",positions.data());
    p.cmd("setMasses",masses.data());
    p.cmd("setBox",box.data());
    p.cmd("setForces",forces.data());
    p.cmd("setVirial",virial.data());
    p.cmd("prepareCalc");
    if(energy) {
      p.cmd("setEnergy",&own);
    }
    p.cmd("performCalc");
    p.cmd("getBias",&bias);
    const double expected=test==3 ? 0 : (energy ? 4 : 1);
    if(std::abs(bias-expected)>1e-12) {
      throw std::runtime_error("Ordinary bias evaluation changed");
    }
    if(!hrex) {
      output<<"ENERGY without HREX: supported\n";
      continue;
    }
    p.cmd("GREX cacheLocalUNow",&own);
    bool rejected=false;
    try {
      p.cmd("GREX cacheLocalUSwap",&foreign);
    } catch(const PLMD::Plumed::Exception& error) {
      if(std::string(error.what()).find("HREX with ENERGY is not supported")==std::string::npos) {
        throw;
      }
      rejected=true;
    }
    if(!energy && rejected) {
      throw std::runtime_error("Coordinate HREX was rejected");
    }
    output<<(test==1 ? "HREX without ENERGY: " : (test==2 ? "HREX with biased ENERGY: " : "HREX with unused ENERGY: "))
          <<(rejected ? "rejected" : "supported")<<'\n';
  }
}
