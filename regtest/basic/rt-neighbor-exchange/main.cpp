#include "plumed/wrapper/Plumed.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Result {
  std::vector<std::array<double,12>> forces;
  std::vector<double> energies;
};

static int configuration(bool accept,int step) {
  int result=0;
  for(int s=0; s<step; ++s)
    if(accept && s%3==1 && s%2==1) {
      result=1-result;
    }
  return result;
}

static Result run(const std::string& name,bool voronoi,int stride,bool trial,bool before,bool accept,bool auxiliary=false) {
  PLMD::Plumed p;
  p.cmd("setNatoms",4);
  p.cmd("setTimestep",.002);
  p.cmd("setLogFile",(name+".log").c_str());
  p.cmd("init");
  const std::string atoms=voronoi?"VORONOI_COORDINATION CENTERS=1-3 ASSIGNED=4 KAPPA=5 REFERENCE=0 SELECT=1"
                          :"COORDINATION GROUPA=1-3 GROUPB=4 R_0=0.3";
  p.cmd("readInputLine",("c: "+atoms+" NOPBC NLIST NL_CUTOFF=1.0 NL_STRIDE="+std::to_string(stride)).c_str());
  p.cmd("readInputLine","BIASVALUE ARG=c");
  std::array<double,12> x{},f{};
  std::array<double,9> virial{};
  const std::array<double,9> box{{100,0,0,0,100,0,0,0,100}};
  const std::array<double,4> masses{{1,1,1,1}};
  auto bind=[&](int step,int conf) {
    x.fill(0);
    x[3]=1.0;
    x[6]=2.0;
    x[9]=conf?1.6:.4;
    p.cmd("setStep",step);
    p.cmd("setPositions",x.data());
    p.cmd("setMasses",masses.data());
    p.cmd("setBox",box.data());
    p.cmd("setForces",f.data());
    p.cmd("setVirial",virial.data());
  };
  auto foreign=[&](int step,int conf) {
    bind(step,1-conf);
    p.cmd("setExchangeStep",1);
    p.cmd("prepareCalc");
    p.cmd("performCalcNoForces");
    p.cmd("setExchangeStep",0);
  };
  Result result;
  for(int step=0; step<12; ++step) {
    const int conf=configuration(accept,step);
    const bool exchange=trial && step%3==1;
    f.fill(0);
    virial.fill(0);
    if(exchange && before) {
      foreign(step,conf);
    }
    bind(step,conf);
    p.cmd("calc");
    double bias=0;
    p.cmd("getBias",&bias);
    result.energies.push_back(bias);
    result.forces.push_back(f);
    if(exchange && !before) {
      foreign(step,conf);
    }
    if(auxiliary)
      for(int i=0; i<3; ++i) {
        bind(step,conf);
        p.cmd("setExchangeStep",1);
        p.cmd("prepareCalc");
        p.cmd("setExchangeStep",0);
        p.cmd("performCalcNoForces");
      }
    if(f!=result.forces.back()) {
      throw std::runtime_error("Read-only calculation changed MD forces");
    }
  }
  return result;
}

static void compare(const Result& reference,const Result& actual,const std::string& name,std::ofstream& output) {
  double error=0;
  for(unsigned step=0; step<reference.forces.size(); ++step) {
    error=std::max(error,std::abs(reference.energies[step]-actual.energies[step]));
    for(unsigned i=0; i<12; ++i) {
      error=std::max(error,std::abs(reference.forces[step][i]-actual.forces[step][i]));
    }
  }
  if(error>1e-10) {
    throw std::runtime_error(name+": invalid neighbor-list energy/forces");
  }
  output<<name<<": pass\n";
}

int main() {
  std::ofstream output("output");
  for(bool voronoi : {
        false,true
      }) {
    const std::string kind=voronoi?"voronoi":"coordination";
    for(bool accept : {
          false,true
        }) {
      const std::string label=kind+(accept?"-accepted":"-rejected");
      const auto reference=run(label+"-reference",voronoi,1,false,false,accept);
      for(bool before : {
            false,true
          }) {
        const std::string name=label+(before?"-before":"-after");
        compare(reference,run(name,voronoi,10,true,before,accept),name,output);
      }
    }
    const auto reference=run(kind+"-fixed",voronoi,1,false,false,false);
    compare(reference,run(kind+"-MCMC",voronoi,10,false,false,false,true),kind+"-MCMC",output);
  }
}
