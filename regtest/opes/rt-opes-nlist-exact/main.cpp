#include "plumed/wrapper/Plumed.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <random>
#include <string>

// OPES with the exact neighbor list must give the same bias and forces as OPES without it
int main() {
  std::ofstream output("output");
  bool all_passed=true;
  //default: adaptive sigma gives kernels of different widths, many of them outside the neighbor list
  //merges: with COMPRESSION_THRESHOLD close to KERNEL_CUTOFF and a small list, a recursive merge
  //can involve a kernel outside the list (name, keywords for both, keywords for the neighbor list)
  const std::array<std::array<std::string,3>,2> settings{{
      {{"default","",""}},{{"merges"," KERNEL_CUTOFF=3.5 COMPRESSION_THRESHOLD=3.4"," NLIST_PARAMETERS=1.1"}}
    }};
  for(const auto& setting : settings)
    for(const std::string mode : {
          "OPES_METAD","OPES_METAD_EXPLORE"
        }) {
    const std::string opes=" ARG=d1,d2 PACE=5 BARRIER=20 TEMP=300"+setting[1];
    PLMD::Plumed p[2];
    int precision=8, natoms=3;
    double timestep=0.001, kbt=2.4943387854;
    for(unsigned n=0; n<2; n++) {
      p[n].cmd("setRealPrecision",&precision);
      p[n].cmd("setNatoms",&natoms);
      p[n].cmd("setTimestep",&timestep);
      p[n].cmd("setKbT",&kbt);
      p[n].cmd("init");
      p[n].cmd("readInputLine","d1: DISTANCE ATOMS=1,2 NOPBC");
      p[n].cmd("readInputLine","d2: DISTANCE ATOMS=1,3 NOPBC");
    }
    const std::string files=mode+"-"+setting[0];
    p[0].cmd("readInputLine",("b: "+mode+opes+" NLIST"+setting[2]+" FILE="+files+"-nlist.kernels").c_str());
    p[1].cmd("readInputLine",("b: "+mode+opes+" FILE="+files+"-ref.kernels").c_str());
    std::mt19937 rng(12345);
    auto uniform=[&]() {
      return rng()/4294967296.0-0.5;
    };
    std::array<double,9> positions{{0,0,0,1,0,0,0,1,0}};
    std::array<double,9> box{{3,0,0,0,3,0,0,0,3}};
    std::array<double,3> masses{{1,1,1}};
    double max_dbias=0, max_dforce=0;
    for(int step=0; step<8000; ++step) {
      for(unsigned j=3; j<9; j++) {
        positions[j]+=0.03*uniform();
      }
      double bias[2];
      std::array<double,9> forces[2];
      for(unsigned n=0; n<2; n++) {
        std::array<double,9> virial{};
        forces[n].fill(0);
        p[n].cmd("setStep",&step);
        p[n].cmd("setPositions",positions.data());
        p[n].cmd("setMasses",masses.data());
        p[n].cmd("setBox",box.data());
        p[n].cmd("setForces",forces[n].data());
        p[n].cmd("setVirial",virial.data());
        p[n].cmd("calc");
        p[n].cmd("getBias",&bias[n]);
      }
      max_dbias=std::max(max_dbias,std::abs(bias[0]-bias[1])/(1+std::abs(bias[1])));
      for(unsigned j=0; j<9; j++) {
        max_dforce=std::max(max_dforce,std::abs(forces[0][j]-forces[1][j])/(1+std::abs(forces[1][j])));
      }
    }
    const bool passed=max_dbias<1e-9 && max_dforce<1e-9;
    all_passed=all_passed && passed;
    output<<mode<<" "<<setting[0]<<": exact neighbor list "<<(passed ? "PASS" : "FAIL")<<"\n";
  }
  return all_passed ? 0 : 1;
}
