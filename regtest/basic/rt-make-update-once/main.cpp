#include "plumed/wrapper/Plumed.h"
#include <array>
#include <fstream>
#include <stdexcept>
#include <string>

void check(long long first_step,const std::string& name) {
  {
    PLMD::Plumed p;
    const int natoms=2;
    std::array<double,6> positions{{0,0,0,1,0,0}},forces{};
    std::array<double,2> masses{{1,1}};
    std::array<double,9> box{{10,0,0,0,10,0,0,0,10}},virial{};
    p.cmd("setNatoms",&natoms);
    p.cmd("setLogFile",(name+".log").c_str());
    p.cmd("init");
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    p.cmd("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1");
    p.cmd("readInputLine",("PRINT ARG=d STRIDE=1 FILE="+name+".colvar").c_str());
    auto evaluate=[&](long long step) {
      forces.fill(0);
      virial.fill(0);
      p.cmd("setStepLongLong",&step);
      p.cmd("setPositions",positions.data());
      p.cmd("setMasses",masses.data());
      p.cmd("setBox",box.data());
      p.cmd("setForces",forces.data());
      p.cmd("setVirial",virial.data());
      p.cmd("prepareCalc");
      p.cmd("performCalcNoUpdate");
    };
    evaluate(first_step);
    p.cmd("update");
    // Amber repeats a complete calculation during startup. Keep every update
    // and output record, warning for each repeated update on that first step.
    evaluate(first_step);
    p.cmd("update");
    evaluate(first_step);
    p.cmd("update");
    // Force-only evaluations must not count as updates or produce warnings.
    evaluate(first_step);
    evaluate(first_step);
    evaluate(first_step+1);
    evaluate(first_step+1);
    p.cmd("update");
    auto expect_error=[&](long long step) {
      bool detected=false;
      try {
        p.cmd("update");
      } catch(const PLMD::Plumed::ExceptionError& error) {
        detected=std::string(error.what()).find("update() called more than once at step "+std::to_string(step))!=std::string::npos;
        if(!detected) {
          throw;
        }
      }
      if(!detected) {
        throw std::runtime_error("duplicate update after the first step did not raise an error");
      }
    };
    expect_error(first_step+1);
    // Returning to the starting step number must not reopen the exception.
    evaluate(first_step);
    p.cmd("update");
    evaluate(first_step);
    expect_error(first_step);
  }
  std::ifstream log(name+".log");
  std::string line;
  unsigned warnings=0;
  while(std::getline(log,line)) {
    if(line.find("WARNING: update() called more than once")!=std::string::npos) {
      if(line.find("at step "+std::to_string(first_step)+",")==std::string::npos) {
        throw std::runtime_error("unexpected duplicate-update warning");
      }
      ++warnings;
    }
  }
  if(warnings!=2) {
    throw std::runtime_error("expected exactly two first-step duplicate-update warnings");
  }
  std::ifstream colvar(name+".colvar");
  unsigned records=0;
  while(std::getline(colvar,line)) {
    if(!line.empty() && line[0]!='#') {
      ++records;
    }
  }
  if(records!=5) {
    throw std::runtime_error("warning must preserve updates; rejected updates must not write records");
  }
}

int main() {
  check(0,"zero");
  check(7,"nonzero");
  check(-1,"negative");
  std::ofstream output("output");
  output << "first-step duplicate updates warn\nsame-step evaluation allowed\nnext step allowed\nlater duplicate updates rejected\nnonzero and negative starts supported\nreturn to initial step remains strict\n";
}
