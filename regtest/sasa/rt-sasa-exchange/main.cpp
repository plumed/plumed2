#include "plumed/wrapper/Plumed.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

static std::string actionKind;
static std::string optionalLibrary;

struct Result {
  std::vector<std::vector<double>> forces;
  std::vector<double> energies;
};

static std::vector<double> structure() {
  std::ifstream pdb("helix.pdb");
  if(!pdb) {
    throw std::runtime_error("Missing helix.pdb");
  }
  std::vector<double> xyz;
  std::string line;
  while(std::getline(pdb,line))
    if(line.compare(0,4,"ATOM")==0) {
      const unsigned atom=std::stoul(line.substr(6,5))-1;
      xyz.resize(std::max(xyz.size(),3*static_cast<std::size_t>(atom+1)));
      for(unsigned d=0; d<3; ++d) {
        xyz[3*atom+d]=std::stod(line.substr(30+8*d,8))*.1;
      }
    }
  return xyz;
}

static Result run(const std::string& name,int stride,bool foreign,bool before,bool accept,bool dense=false) {
  PLMD::Plumed p;
  const auto original=structure();
  const int natoms=original.size()/3;
  p.cmd("setNatoms",natoms);
  p.cmd("setTimestep",.002);
  p.cmd("setLogFile",(name+".log").c_str());
  p.cmd("init");
  if(!optionalLibrary.empty()) {
    p.cmd("readInputLine",("LOAD FILE="+optionalLibrary).c_str());
  }
  p.cmd("readInputLine","MOLINFO MOLTYPE=protein STRUCTURE=helix.pdb");
  p.cmd("readInputLine",("s: "+actionKind+" ATOMS=7-126 TYPE=TOTAL NOPBC NL_STRIDE="+std::to_string(stride)).c_str());
  p.cmd("readInputLine","BIASVALUE ARG=s");
  std::vector<double> x(original.size()),f(original.size()),masses(natoms,1.0);
  const std::array<double,9> box{{100,0,0,0,100,0,0,0,100}};
  std::array<double,9> virial{};
  auto bind=[&](int step,int conf) {
    for(unsigned i=0; i<x.size(); ++i) {
      x[i]=original[i]*(conf?4.0:1.0);
    }
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
  Result result;
  for(int step=0; step<5; ++step) {
    const int conf=accept && step>=2 ? 1 : 0;
    std::fill(f.begin(),f.end(),0.0);
    virial.fill(0.0);
    const bool exchange=foreign && (dense || step==1);
    if(exchange && before) {
      trial(step,conf);
    }
    bind(step,conf);
    p.cmd("calc");
    double bias=0;
    p.cmd("getBias",&bias);
    result.energies.push_back(bias);
    result.forces.push_back(f);
    if(exchange && !before) {
      trial(step,conf);
      if(f!=result.forces.back()) {
        throw std::runtime_error("Foreign evaluation changed MD forces");
      }
    }
  }
  return result;
}

int main(int argc,char** argv) {
  if(argc>1) {
    optionalLibrary=argv[1];
  }
  std::ofstream output("output");
  for(const auto& kind : {
        "SASA_LCPO","SASA_HASEL"
      }) {
    actionKind=kind;
    for(bool accept : {
          false,true
        }) {
      const std::string label=actionKind+"-"+(accept?"accepted":"rejected");
      const auto reference=run(label+"-reference",1,false,false,accept);
      for(bool dense : {
            false,true
          })
        for(bool before : {
              false,true
            }) {
          const std::string name=label+(dense?"-dense":"-sparse")+(before?"-before":"-after");
          const auto actual=run(name,40,true,before,accept,dense);
          double error=0.0;
          for(unsigned step=0; step<reference.forces.size(); ++step) {
            error=std::max(error,std::abs(reference.energies[step]-actual.energies[step]));
            for(unsigned atom=0; atom<reference.forces[step].size(); ++atom) {
              error=std::max(error,std::abs(reference.forces[step][atom]-actual.forces[step][atom]));
            }
          }
          if(error>1e-10) {
            throw std::runtime_error(name+": foreign neighbor list changed local energies/forces; error="+std::to_string(error));
          }
          output<<name<<": pass\n";
        }
    }
  }
}
