#ifdef __PLUMED_LIBCXX11
// In order to correctly catch the thrown C++11 exceptions,
// we notify the Plumed wrapper that those exceptions are recognized by the compiler.
#define __PLUMED_WRAPPER_LIBCXX11 1
#endif
#include "plumed/wrapper/Plumed.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// The final checkpoint (CHECKPOINT, plumed_final.cpt) written when PLUMED is
// finalized by an MD engine that never sends the checkpointing signal and does
// not call runFinalJobs, and the cases in which it is not written.

using namespace PLMD;

namespace {

const int natoms=4;

std::string readFile(const std::string& name) {
  std::ifstream f(name);
  std::stringstream ss;
  ss<<f.rdbuf();
  return ss.str();
}

bool exists(const std::string& name) {
  return std::ifstream(name).good();
}

// the lines of the log of a run that contain this text (without the PLUMED: prefix)
std::string grepLog(const std::string& log,const std::string& text) {
  std::ifstream f(log);
  std::string line,out;
  while(std::getline(f,line)) {
    if(line.find(text)!=std::string::npos) {
      out+="  "+line.substr(line.find(':')+2)+"\n";
    }
  }
  return out;
}

struct Engine {
  Plumed p;
  std::vector<double> pos,forces,masses;
  double box[9]= {10,0,0,0,10,0,0,0,10};
  double virial[9]= {0};
  Engine(const std::string& name,bool restart,long long int restartStep=-1):
    pos(3*natoms),forces(3*natoms),masses(natoms,1.0) {
    int n=natoms;
    p.cmd("setNatoms",&n);
    p.cmd("setMDEngine","test");
    p.cmd("setTimestep",0.1);
    p.cmd("setLogFile",(name+".log").c_str());
    if(restart) {
      p.cmd("setRestart",1);
      p.cmd("setRestartStep",restartStep);
    }
    p.cmd("init");
    p.cmd("readInputLine",("CHECKPOINT FILE="+name+".cpt").c_str());
    p.cmd("readInputLine","d1: DISTANCE ATOMS=1,2");
    p.cmd("readInputLine","d2: DISTANCE ATOMS=3,4");
    // a state that goes from one step to the next (the work)
    p.cmd("readInputLine","mr: MOVINGRESTRAINT ARG=d1,d2 STEP0=0 AT0=1.0,1.0 KAPPA0=100,100 STEP1=60 AT1=2.0,0.5 KAPPA1=200,100");
    p.cmd("readInputLine",("PRINT ARG=* FILE=COLVAR_"+name+" FMT=%8.4f").c_str());
  }
  void prepare(long long int step) {
    for(int i=0; i<natoms; i++) {
      pos[3*i]=i+0.5*std::sin(0.1*step+i);
      pos[3*i+1]=0.3*std::cos(0.07*step*(i+1));
      pos[3*i+2]=0.2*i*step/60.0;
    }
    p.cmd("setStepLongLong",step);
    p.cmd("setPositions",pos.data(),3*natoms);
    p.cmd("setForces",forces.data(),3*natoms);
    p.cmd("setMasses",masses.data(),natoms);
    p.cmd("setBox",box,9);
    p.cmd("setVirial",virial,9);
  }
  void step(long long int step,bool signal=false) {
    prepare(step);
    if(signal) {
      int one=1;
      p.cmd("doCheckPoint",&one);
    }
    p.cmd("calc");
  }
};

}

int main() {
  std::ofstream ofs("output");

  // reference: steps 0 to 59 without interruption
  {
    Engine e("ref",false);
    for(int s=0; s<60; s++) {
      e.step(s);
    }
  }

  // A) steps 0 to 24, no signal, no runFinalJobs: the final checkpoint (for step 25) is written
  // when PLUMED is finalized; the restart from it at step 25 reproduces the reference
  {
    Engine e("a",false);
    for(int s=0; s<25; s++) {
      e.step(s);
    }
  }
  ofs<<"A: a_final.cpt "<<(exists("a_final.cpt")?"written":"not written")<<"\n"<<grepLog("a.log","inal checkpoint");
  {
    Engine e("a",true,25);
    for(int s=25; s<60; s++) {
      e.step(s);
    }
  }
  ofs<<grepLog("a.log","reading the checkpoint")<<grepLog("a.log","restored");
  ofs<<"A: COLVAR "<<(readFile("COLVAR_a")==readFile("COLVAR_ref")?"identical":"different")<<" to the reference\n";

  // B) runFinalJobs, then finalization: written only once
  {
    Engine e("b",false);
    for(int s=0; s<25; s++) {
      e.step(s);
    }
    e.p.cmd("runFinalJobs");
  }
  ofs<<"B:\n"<<grepLog("b.log","b_final.cpt");

  // C) the checkpointing signal at step 10: no final checkpoint
  {
    Engine e("c",false);
    for(int s=0; s<25; s++) {
      e.step(s,s==10);
    }
  }
  ofs<<"C: c_final.cpt "<<(exists("c_final.cpt")?"written":"not written")<<"\n";

  // D) the last step is not completed (calculation without update): not written
  {
    Engine e("d",false);
    for(int s=0; s<25; s++) {
      e.step(s);
    }
    e.prepare(25);
    e.p.cmd("prepareCalc");
    e.p.cmd("performCalcNoUpdate");
  }
  ofs<<"D: d_final.cpt "<<(exists("d_final.cpt")?"written":"not written")<<"\n"<<grepLog("d.log","inal checkpoint");

  // E) PLUMED finalized while an exception is propagating: not written
  try {
    Engine e("e",false);
    for(int s=0; s<25; s++) {
      e.step(s);
    }
    throw std::runtime_error("the MD engine failed");
  } catch(const std::runtime_error&) {
  }
  ofs<<"E: e_final.cpt "<<(exists("e_final.cpt")?"written":"not written")<<"\n";
  return 0;
}
