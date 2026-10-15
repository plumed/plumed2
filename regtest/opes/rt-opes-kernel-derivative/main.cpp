#include "plumed/wrapper/Plumed.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <string>

int main() {
  std::ofstream output("output");
  for(const std::string mode : {
        "OPES_METAD","OPES_METAD_EXPLORE"
      }) {
    PLMD::Plumed p;
    int precision=8, natoms=2;
    double timestep=0.001, kbt=2.4943387854;
    p.cmd("setRealPrecision",&precision);
    p.cmd("setNatoms",&natoms);
    p.cmd("setTimestep",&timestep);
    p.cmd("setKbT",&kbt);
    p.cmd("init");
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    p.cmd("readInputLine",("b: "+mode+
                           " ARG=d PACE=1 BARRIER=10 TEMP=300 SIGMA=0.1 KERNEL_CUTOFF=3 FILE="+mode+".kernels").c_str());
    std::array<double,6> positions{{0,0,0,0.5,0.2,0.1}}, forces{};
    std::array<double,9> box{{3,0,0,0,3,0,0,0,3}}, virial{};
    std::array<double,2> masses{{1,1}};
    int step=0;
    auto evaluate = [&](bool update) {
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
      return bias;
    };
    for(step=0; step<4; ++step) {
      positions[3]=0.5+0.02*step;
      evaluate(true);
    }
    positions[3]=0.53;
    evaluate(false);
    const auto analytic=forces;
    double maximum=0;
    const double h=1e-6;
    for(unsigned i=0; i<positions.size(); ++i) {
      const double x=positions[i];
      positions[i]=x+h;
      const double plus=evaluate(false);
      positions[i]=x-h;
      const double minus=evaluate(false);
      positions[i]=x;
      const double numerical=-(plus-minus)/(2*h);
      maximum=std::max(maximum,std::abs(analytic[i]));
      if(!std::isfinite(numerical) ||
          std::abs(analytic[i]-numerical)>1e-7*(1+std::abs(numerical))) {
        output<<mode<<" derivative FAIL "<<i<<" "<<analytic[i]<<" "<<numerical<<"\n";
        return 1;
      }
    }
    if(maximum<1e-3) {
      return 2;
    }
    output<<mode<<" derivative PASS\n";
  }
}
