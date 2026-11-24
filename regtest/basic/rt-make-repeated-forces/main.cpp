#define __PLUMED_WRAPPER_IMPLEMENTATION 1
#define __PLUMED_WRAPPER_LINK_RUNTIME 1
#include "plumed/wrapper/Plumed.h"
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

struct Engine {
  PLMD::Plumed p;
  std::array<double,6> x{},f{};
  std::array<double,2> m{{1,1}},q{};
  std::array<double,9> v{},box{{10,0,0,0,10,0,0,0,10}};
  explicit Engine(const std::string& input) {
    p.cmd("setNatoms",2);
    p.cmd("setTimestep",0.001);
    p.cmd("setLogFile","engine.log");
    p.cmd("init");
    p.cmd("readInputLines",input.c_str());
  }
  double evaluate(int step,double distance,const char* command) {
    x[3]=distance;
    f.fill(0);
    v.fill(0);
    p.cmd("setStep",step);
    p.cmd("setPositions",x.data());
    p.cmd("setMasses",m.data());
    p.cmd("setCharges",q.data());
    p.cmd("setBox",box.data());
    p.cmd("setForces",f.data());
    p.cmd("setVirial",v.data());
    p.cmd("prepareCalc");
    p.cmd(command);
    double bias=0;
    p.cmd("getBias",&bias);
    return bias;
  }
};
void near(double actual,double expected) {
  if(!std::isfinite(actual) || std::abs(actual-expected)>1e-12) {
    throw std::runtime_error("repeated evaluation changed bias/state or lost forces: actual="+std::to_string(actual)+" expected="+std::to_string(expected));
  }
}
int main() {
  const std::string guest="d: DISTANCE ATOMS=1,2 NOPBC\na: ABMD ARG=d TO=0 KAPPA=2\nr: RESTRAINT ARG=d AT=0 KAPPA=2\n";
  std::ofstream("guest.dat")<<guest;
  const char* kernel=std::getenv("PLUMED_KERNEL");
  if(!kernel) {
    throw std::runtime_error("PLUMED_KERNEL is required for a loaded-guest comparison");
  }
  for(int mode=0; mode<3; ++mode) {
    for(bool split : {
          false,true
        }) {
      std::string input=mode==0 ? guest : "p: PLUMED FILE=guest.dat";
      if(mode==2) {
        input+=" KERNEL="+std::string(kernel);
      }
      std::cerr<<"mode="<<mode<<" split="<<split<<"\n";
      Engine engine(input);
      near(engine.evaluate(0,2,split ? "performCalcNoUpdate" : "performCalc"),4);
      near(engine.f[0],4);
      if(split) {
        engine.p.cmd("update");
      }
      engine.p.cmd("setEvaluationMode",1);
      near(engine.evaluate(0,1,"performCalcNoUpdate"),1);
      engine.p.cmd("setEvaluationMode",0);
      near(engine.f[0],2);
      near(engine.f[3],-2);
      near(engine.evaluate(1,2,"performCalc"),4);
      near(engine.f[0],4);
    }
  }
  const std::string strided="d: DISTANCE ATOMS=1,2 NOPBC\nr: RESTRAINT ARG=d AT=0 KAPPA=2 STRIDE=2\nPRINT ARG=d FILE=STRIDE_COLVAR\n";
  std::ofstream("strided.dat")<<strided;
  for(int mode=0; mode<3; ++mode) {
    std::string input=mode==0 ? strided : "p: PLUMED FILE=strided.dat";
    if(mode==2) {
      input+=" KERNEL="+std::string(kernel);
    }
    Engine engine(input);
    near(engine.evaluate(0,2,"performCalc"),4);
    near(engine.f[0],8);
    near(engine.evaluate(1,2,"performCalc"),0);
    near(engine.f[0],0);
    engine.p.cmd("setEvaluationMode",1);
    near(engine.evaluate(1,2,"performCalcNoUpdate"),0);
    near(engine.f[0],0);
    engine.p.cmd("setEvaluationMode",2);
    near(engine.evaluate(1,2,"performCalcNoForces"),4);
    near(engine.f[0],0);
    engine.p.cmd("setEvaluationMode",0);
    near(engine.evaluate(2,2,"performCalc"),4);
    near(engine.f[0],8);
  }
  {
    Engine metad("d: DISTANCE ATOMS=1,2 NOPBC\nm: METAD ARG=d SIGMA=1 HEIGHT=1 PACE=1 FILE=HILLS\n");
    near(metad.evaluate(0,2,"performCalc"),0);
    near(metad.evaluate(1,2,"performCalc"),0);
    near(metad.evaluate(1,2,"performCalc"),1);
    near(metad.evaluate(1,2,"performCalc"),2);
    near(metad.evaluate(2,2,"performCalc"),3);
  }
  std::ifstream hills("HILLS");
  std::string line;
  int deposits=0;
  while(std::getline(hills,line)) {
    if(!line.empty() && line[0]!='#') {
      ++deposits;
    }
  }
  if(deposits!=4) {
    throw std::runtime_error("normal-mode repeated full calculations changed legacy semantics");
  }
  std::ofstream("output")<<"Explicit auxiliary force queries preserve ABMD and nested forces; normal commands retain semantics: OK\n";
}
