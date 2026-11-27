#define __PLUMED_WRAPPER_IMPLEMENTATION 1
#define __PLUMED_WRAPPER_LINK_RUNTIME 1
#include "plumed/wrapper/Plumed.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<double> run(int mode) {
  const std::string colvar="COLVAR"+std::to_string(mode);
  std::ofstream("acc.in")<<"#! FIELDS time m.acc\n0.009 2\n";
  {
    PLMD::Plumed p;
    std::array<double,6> x{{0,0,0,2,0,0}},f{};
    std::array<double,2> masses{{1,1}};
    std::array<double,9> box{{10,0,0,0,10,0,0,0,10}},virial{};
    p.cmd("setNatoms",2);
    p.cmd("setTimestep",0.001);
    p.cmd("setLogFile","plumed.log");
    p.cmd("init");
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    p.cmd("readInputLine","m: METAD ARG=d SIGMA=1 HEIGHT=1 PACE=1000 TEMP=300 ACCELERATION ACCELERATION_RFILE=acc.in FILE=HILLS"+std::to_string(mode));
    p.cmd("readInputLine","PRINT ARG=m.acc FILE="+colvar+" FMT=%.17g");
    auto bind=[&](int step) {
      f.fill(0);
      virial.fill(0);
      p.cmd("setStep",step);
      p.cmd("setPositions",x.data());
      p.cmd("setMasses",masses.data());
      p.cmd("setBox",box.data());
      p.cmd("setForces",f.data());
      p.cmd("setVirial",virial.data());
    };
    bind(10);
    p.cmd("calc");
    if(mode) {
      bind(10);
      if(mode==2) {
        p.cmd("setExchangeStep",1);
      }
      if(mode==3 || mode==4) {
        p.cmd("setEvaluationMode",mode==3 ? 1 : 2);
      }
      p.cmd("prepareCalc");
      p.cmd(mode==3 ? "performCalcNoUpdate" : "performCalcNoForces");
      if(mode==2) {
        p.cmd("setExchangeStep",0);
      }
      if(mode==3 || mode==4) {
        p.cmd("setEvaluationMode",0);
      }
    }
    bind(11);
    p.cmd("calc");
  }
  std::ifstream file(colvar);
  std::string line;
  std::vector<double> values;
  while(std::getline(file,line)) {
    if(line.empty() || line[0]=='#') {
      continue;
    }
    std::istringstream row(line);
    double time,value;
    row>>time>>value;
    values.push_back(value);
  }
  if(values.size()!=2) {
    throw std::runtime_error("missing acceleration samples");
  }
  return values;
}
int main() {
  PLMD::Plumed probe;
  int api=0;
  probe.cmd("getApiVersion",&api);
  const auto control=run(0);
  for(int mode : {
        1,2,3,4
      }) {
    const auto trial=run(mode);
    std::cerr<<"mode="<<mode<<" next="<<trial[1]<<" control="<<control[1]<<"\n";
    for(unsigned i=0; i<control.size(); ++i) {
      if(!std::isfinite(trial[i]) || std::abs(trial[i]-control[i])>1e-12) {
        throw std::runtime_error("auxiliary evaluation changed restarted acceleration");
      }
    }
  }
  std::ofstream("output")<<"Restarted METAD acceleration is unchanged by auxiliary evaluations: OK\n";
}
