#include "plumed/wrapper/Plumed.h"
#include "plumed/config/Config.h"
#include <array>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

int main() {
  int major=0,minor=0;
  char dot=0;
  std::istringstream version(PLMD::config::getVersion());
  version >> major >> dot >> minor;
  const bool expect_error=major>2 || (major==2 && minor>=11);
  bool detected=false;
  {
    PLMD::Plumed p;
    const int natoms=2;
    std::array<double,6> positions{{0,0,0,1,0,0}},forces{};
    std::array<double,2> masses{{1,1}};
    std::array<double,9> box{{10,0,0,0,10,0,0,0,10}},virial{};
    p.cmd("setNatoms",&natoms);
    p.cmd("setLogFile","update.log");
    p.cmd("init");
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    p.cmd("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1");
    auto evaluate=[&](int step) {
      forces.fill(0);
      virial.fill(0);
      p.cmd("setStep",&step);
      p.cmd("setPositions",positions.data());
      p.cmd("setMasses",masses.data());
      p.cmd("setBox",box.data());
      p.cmd("setForces",forces.data());
      p.cmd("setVirial",virial.data());
      p.cmd("prepareCalc");
      p.cmd("performCalcNoUpdate");
    };
    evaluate(0);
    p.cmd("update");
    evaluate(0);
    evaluate(0);
    evaluate(1);
    p.cmd("update");
    try {
      p.cmd("update");
    } catch(const PLMD::Plumed::ExceptionError& error) {
      detected=std::string(error.what()).find("update() called more than once at step 1")!=std::string::npos;
      if(!expect_error || !detected) {
        throw;
      }
    }
  }
  std::ifstream log("update.log");
  std::string line;
  unsigned warnings=0;
  while(std::getline(log,line)) {
    if(line.find("WARNING: update() called more than once")!=std::string::npos) {
      if(line.find("at step 1")==std::string::npos || line.find("error in PLUMED 2.11")==std::string::npos) {
        throw std::runtime_error("unexpected duplicate-update warning");
      }
      ++warnings;
    }
  }
  if(expect_error ? (!detected || warnings!=0) : (detected || warnings!=1)) {
    throw std::runtime_error("duplicate update did not follow the release policy");
  }
  std::ofstream output("output");
  output << "same-step evaluation allowed\nnext step allowed\nduplicate update detected\n";
}
