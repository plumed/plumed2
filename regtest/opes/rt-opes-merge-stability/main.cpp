#include "plumed/wrapper/Plumed.h"
#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

struct Result {
  std::vector<double> samples;
  bool finite=true;
};

Result evaluate(const std::string& mode, bool shifted, bool moving) {
  Result result;
  PLMD::Plumed p;
  int precision=8, natoms=2;
  double timestep=0.001, kbt=2.4943387854;
  p.cmd("setRealPrecision",&precision);
  p.cmd("setNatoms",&natoms);
  p.cmd("setTimestep",&timestep);
  p.cmd("setKbT",&kbt);
  p.cmd("init");
  p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
  p.cmd("readInputLine",shifted ?
        "c: CUSTOM ARG=d FUNC=x+100000000 PERIODIC=NO" :
        "c: CUSTOM ARG=d FUNC=x PERIODIC=NO");
  const std::string name=mode+(shifted ? "-shifted" : "-unshifted")+
                         (moving ? "-moving" : "-constant");
  p.cmd("readInputLine",("b: "+mode+
        " ARG=c PACE=1 BARRIER=10 TEMP=300 SIGMA=0.125 FIXED_SIGMA "
        "COMPRESSION_THRESHOLD=1 KERNEL_CUTOFF=3 FILE="+name+".kernels").c_str());
  std::array<double,6> positions{{0,0,0,1,0,0}}, forces{};
  std::array<double,9> box{{3,0,0,0,3,0,0,0,3}}, virial{};
  std::array<double,2> masses{{1,1}};
  int step=0;
  auto sample=[&](bool update) {
    forces.fill(0);
    virial.fill(0);
    p.cmd("setStep",&step);
    p.cmd("setPositions",positions.data());
    p.cmd("setMasses",masses.data());
    p.cmd("setBox",box.data());
    p.cmd("setForces",forces.data());
    p.cmd("setVirial",virial.data());
    if(update) {
      p.cmd("calc");
    } else {
      p.cmd("prepareCalc");
      p.cmd("performCalcNoUpdate");
    }
    double bias=0;
    p.cmd("getBias",&bias);
    result.samples.push_back(bias);
    result.samples.insert(result.samples.end(),forces.begin(),forces.end());
    result.samples.insert(result.samples.end(),virial.begin(),virial.end());
  };
  for(step=0; step<4; ++step) {
    positions[3]=1+(moving ? (step%2)*0.03125 : 0);
    sample(true);
  }
  positions[3]=1.0625;
  sample(false);
  for(double value : result.samples) {
    result.finite=result.finite && std::isfinite(value);
  }
  return result;
}

int main() {
  std::ofstream output("output"), details("details");
  bool passed=true;
  for(const std::string mode : {"OPES_METAD","OPES_METAD_EXPLORE"}) {
    for(bool moving : {false,true}) {
      const Result reference=evaluate(mode,false,moving);
      const Result shifted=evaluate(mode,true,moving);
      details << mode << " moving=" << moving
              << " reference_finite=" << reference.finite
              << " shifted_finite=" << shifted.finite << "\n";
      for(unsigned i=0; i<reference.samples.size(); ++i) {
        details << i << " " << reference.samples[i] << " "
                << shifted.samples[i] << "\n";
      }
      bool same=reference.finite && shifted.finite;
      for(unsigned i=0; i<reference.samples.size(); ++i) {
        same=same && std::abs(reference.samples[i]-shifted.samples[i])<
             2e-6*(1+std::abs(reference.samples[i]));
      }
      output << mode << (moving ? " moving" : " constant")
             << " translation-invariance " << (same ? "PASS" : "FAIL") << "\n";
      passed=passed && same;
    }
  }
  return passed ? 0 : 1;
}
