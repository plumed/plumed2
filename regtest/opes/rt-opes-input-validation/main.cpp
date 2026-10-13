#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <string>

// Invalid widths and deposition strides must fail at input parsing, before
// their first division or modulo in a simulation.
int main() {
  std::ofstream output("output");
  bool passed=true;
  for(const std::string mode : {
        "OPES_METAD", "OPES_METAD_EXPLORE"
      }) {
    for(const std::string parameters : {
          "PACE=0 SIGMA=0.25",
          "PACE=-1 SIGMA=0.25",
          "PACE=1 SIGMA=0",
          "PACE=1 SIGMA=-0.25",
          "PACE=1 SIGMA=0.25 SIGMA_MIN=0",
          "PACE=1 SIGMA=0.25 SIGMA_MIN=-0.1"
        }) {
      bool rejected=false;
      try {
        PLMD::Plumed p;
        const int natoms=2;
        const double timestep=0.001, kbt=2.49433863;
        p.cmd("setNatoms",&natoms);
        p.cmd("setTimestep",&timestep);
        p.cmd("setKbT",&kbt);
        p.cmd("init");
        p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
        const std::string input="b: "+mode+" ARG=d TEMP=300 BARRIER=10 "+parameters;
        p.cmd("readInputLine",input.c_str());
      } catch(const std::exception& e) {
        rejected=std::string(e.what()).find("greater than zero")!=std::string::npos;
      }
      output << mode << " " << parameters << " " << (rejected ? "PASS" : "FAIL") << "\n";
      passed=passed && rejected;
    }
  }
  return passed ? 0 : 1;
}
