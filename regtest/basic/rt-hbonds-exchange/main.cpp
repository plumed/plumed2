#include "plumed/wrapper/Plumed.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
struct Result {
  std::vector<double> biases;
  std::vector<std::array<double,18>> forces;
};
static Result run(const std::string& name,int stride,bool trials,bool before,bool accept,bool dense) {
  PLMD::Plumed p;
  p.cmd("setNatoms",6);
  p.cmd("setTimestep",.002);
  p.cmd("setLogFile",(name+".log").c_str());
  p.cmd("init");
  p.cmd("readInputLine",("h: HBONDS DONORS=1,4 ACCEPTORS=3,6 HYDROGENS=2,5 NLIST NL_CUTOFF=0.5 NL_STRIDE="+std::to_string(stride)+" NOPBC").c_str());
  p.cmd("readInputLine","BIASVALUE ARG=h");
  std::array<double,18> x{},f{};
  std::array<double,6> masses{{1,1,1,1,1,1}};
  std::array<double,9> box{{10,0,0,0,10,0,0,0,10}},virial{};
  auto bind=[&](int step,int conf) {
    x.fill(0);
    x[3]=.1;
    x[6]=conf?1:.25;
    x[9]=2;
    x[12]=2.1;
    x[15]=conf?3:2.25;
    p.cmd("setStep",step);
    p.cmd("setPositions",x.data());
    p.cmd("setMasses",masses.data());
    p.cmd("setBox",box.data());
    p.cmd("setForces",f.data());
    p.cmd("setVirial",virial.data());
  };
  auto trial=[&](int step,int conf) {
    bind(step,1-conf);
    p.cmd("setExchangeStep",1);
    p.cmd("prepareCalc");
    p.cmd("performCalcNoForces");
    p.cmd("setExchangeStep",0);
  };
  Result r;
  for(int step=0; step<5; ++step) {
    int conf=accept&&step>=2;
    bool attempt=trials&&step>0&&(dense||step==1);
    f.fill(0);
    virial.fill(0);
    if(attempt&&before) {
      trial(step,conf);
    }
    bind(step,conf);
    p.cmd("calc");
    double bias=0;
    p.cmd("getBias",&bias);
    r.biases.push_back(bias);
    r.forces.push_back(f);
    if(attempt&&!before) {
      trial(step,conf);
    }
    if(f!=r.forces.back()) {
      throw std::runtime_error("Foreign calculation changed MD forces");
    }
  }
  return r;
}
int main() {
  std::ofstream out("output");
  for(bool accept: {
        false,true
      })
    for(bool dense: {
          false,true
        })
      for(bool before: {
            false,true
          }) {
        const std::string name=std::string(accept?"accepted":"rejected")+(dense?"-dense":"-sparse")+(before?"-before":"-after");
        const auto reference=run(name+"-reference",1,false,false,accept,false),actual=run(name,40,true,before,accept,dense);
        double error=0;
        for(unsigned s=0; s<5; ++s) {
          error=std::max(error,std::abs(reference.biases[s]-actual.biases[s]));
          for(unsigned i=0; i<18; ++i) {
            error=std::max(error,std::abs(reference.forces[s][i]-actual.forces[s][i]));
          }
        }
        if(error>=1e-10) {
          throw std::runtime_error(name+": local bias/forces changed after trial; error="+std::to_string(error));
        }
        out<<name<<": pass\n";
      }
  return 0;
}
