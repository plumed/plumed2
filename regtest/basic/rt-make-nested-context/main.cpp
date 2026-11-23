#include "plumed/bias/Bias.h"
#include "plumed/core/ActionRegister.h"
#include "plumed/wrapper/Plumed.h"
#include <array>
#include <fstream>
#include <stdexcept>

namespace PLMD {
class TrialContextProbe : public bias::Bias {
public:
  static void registerKeywords(Keywords& keys) {
    bias::Bias::registerKeywords(keys);
  }
  explicit TrialContextProbe(const ActionOptions& options):Action(options),Bias(options) {}
  void calculate() override {
    setBias(getExchangeStep() ? getArgument(0) : 0.);
    setOutputForce(0,0.);
  }
};
PLUMED_REGISTER_ACTION(TrialContextProbe,"TRIAL_CONTEXT_PROBE")
}

int main() {
  std::ofstream("guest.dat")<<"d: DISTANCE ATOMS=1,2 NOPBC\nc: TRIAL_CONTEXT_PROBE ARG=d\n";
  PLMD::Plumed p;
  p.cmd("setNatoms",2);
  p.cmd("setTimestep",.001);
  p.cmd("setLogFile","plumed.log");
  p.cmd("init");
  p.cmd("readInputLine","g: PLUMED FILE=guest.dat");
  std::array<double,6> x{{0,0,0,2,0,0}},f{};
  std::array<double,2> mass{{1,1}},charge{};
  std::array<double,9> box{{10,0,0,0,10,0,0,0,10}},virial{};
  auto bind=[&]() {
    p.cmd("setStep",0);
    p.cmd("setPositions",x.data());
    p.cmd("setMasses",mass.data());
    p.cmd("setCharges",charge.data());
    p.cmd("setForces",f.data());
    p.cmd("setBox",box.data());
    p.cmd("setVirial",virial.data());
  };
  bind();
  p.cmd("calc");
  bind();
  p.cmd("setExchangeStep",1);
  p.cmd("prepareCalc");
  // Matches the GROMACS lambda-trial sequence: activate dependencies, then
  // evaluate without declaring a foreign-coordinate exchange.
  p.cmd("setExchangeStep",0);
  p.cmd("performCalcNoForces");
  double bias=0.;
  p.cmd("getBias",&bias);
  if(bias!=0.) {
    throw std::runtime_error("nested guest retained preparation's exchange flag");
  }
  std::ofstream("output")<<"Nested lambda trial: current exchange context propagated OK\n";
}
