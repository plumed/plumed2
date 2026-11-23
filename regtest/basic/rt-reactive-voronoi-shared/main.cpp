#include "plumed/wrapper/Plumed.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

struct System {
  PLMD::Plumed p;
  std::vector<double> x{0.1,0.2,0.1, 1.5,0.3,0.2, 0.4,0.3,0.2, 1.1,0.5,0.2};
  std::vector<double> f=std::vector<double>(12,0), m=std::vector<double>(4,1);
  std::array<double,9> box{{3,0,0,0,3,0,0,0,3}}, vir{};
  double energy=0;
  System(bool share, const std::string& mode, const std::string& backend="") {
    int precision=8,n=4;
    double dt=0.0005;
    p.cmd("setRealPrecision",&precision);
    p.cmd("setNatoms",&n);
    p.cmd("setTimestep",&dt);
    p.cmd("setLogFile","/dev/null");
    p.cmd("init");
    const std::string extra=(share ? "" : " NO_SHARE")+mode+backend;
    auto read=[&](const std::string& s) {
      p.cmd("readInputLine",s.c_str());
    };
    read("a: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=0.8,1.2 POWER=2"+extra);
    read("b: VORONOI_DISTANCE CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=1 GROUP1=1 GROUP2=2"+extra);
    read("c: VORONOI_POSITION CENTERS=1,2 ASSIGNED=3,4 KAPPA=8 REFERENCE=0.9,1.1 AXIS=Y ORIGIN=0.3"+extra);
    read("d: VORONOI_COORDINATION CENTERS=2,1 ASSIGNED=4,3 KAPPA=5 REFERENCE=1.2,0.8 SELECT=1 COEFFICIENTS=1.3 POWER=2"+extra);
    read("bias: RESTRAINT ARG=a,b,c,d AT=0.3,-0.1,0.2,0.05 KAPPA=3,2,1,4");
  }
  void eval(int step=7) {
    std::fill(f.begin(),f.end(),0);
    vir.fill(0);
    p.cmd("setStep",&step);
    p.cmd("setPositions",x.data());
    p.cmd("setMasses",m.data());
    p.cmd("setBox",box.data());
    p.cmd("setForces",f.data());
    p.cmd("setVirial",vir.data());
    p.cmd("calc");
    p.cmd("getBias",&energy);
  }
};

int main(int argc,char** argv) {
  std::ofstream out("output");
  double maxDifference=0,maxFiniteDifference=0;
  const std::string backend=argc>1 ? argv[1] : "";
  for(const std::string mode: {
        "", " NOPBC", " NLIST NL_CUTOFF=2 NL_STRIDE=1",
        " NLIST NL_CUTOFF=2 NL_STRIDE=3 NL_SKIN=0.4"
      }) {
    if(!backend.empty() && mode.find("NLIST")!=std::string::npos) {
      continue;
    }
    System shared(true,mode,backend), independent(false,mode);
    for(int k=0; k<12; ++k) {
      // Repeated step, changed coordinates and box, then skipped steps.
      shared.x[6]+=0.004;
      shared.box[0]+=0.003;
      independent.x=shared.x;
      independent.box=shared.box;
      const int step=k<8 ? 7 : 7+3*k;
      shared.eval(step);
      independent.eval(step);
      maxDifference=std::max(maxDifference,std::abs(shared.energy-independent.energy));
      for(unsigned i=0; i<12; ++i) {
        maxDifference=std::max(maxDifference,std::abs(shared.f[i]-independent.f[i]));
      }
      for(unsigned i=0; i<9; ++i) {
        maxDifference=std::max(maxDifference,std::abs(shared.vir[i]-independent.vir[i]));
      }
    }
    shared.eval();
    const auto force=shared.f;
    const double h=1e-6;
    for(unsigned i=0; i<12; ++i) {
      const double x=shared.x[i];
      shared.x[i]=x+h;
      shared.eval();
      const double plus=shared.energy;
      shared.x[i]=x-h;
      shared.eval();
      const double minus=shared.energy;
      shared.x[i]=x;
      maxFiniteDifference=std::max(maxFiniteDifference,std::abs(force[i]+(plus-minus)/(2*h)));
    }
  }
  std::cout<<"max_shared_force_energy_virial_error "<<maxDifference
           <<"\nmax_bias_force_finite_difference_error "<<maxFiniteDifference<<"\n";
  const bool pass=maxDifference<2e-11 && maxFiniteDifference<2e-7;
  out<<(pass ? "PASS" : "FAIL")<<" shared multi-bias same-step coordinate/box invalidation and finite differences\n";
  return pass ? 0 : 1;
}
