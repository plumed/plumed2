#include "plumed/wrapper/Plumed.h"
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

struct Evaluation {
  std::string error;
  std::vector<double> forces;
};

Evaluation evaluate(const std::vector<std::string>& input,
                    std::vector<double> positions) {
  Evaluation result;
  result.forces.assign(positions.size(),0.0);
  try {
    PLMD::Plumed p;
    int precision=8, natoms=positions.size()/3;
    double timestep=0.001;
    p.cmd("setRealPrecision",&precision);
    p.cmd("setNatoms",&natoms);
    p.cmd("setTimestep",&timestep);
    p.cmd("init");
    for(const auto& line : input) {
      p.cmd("readInputLine",line.c_str());
    }
    std::vector<double> masses(natoms,1.0);
    std::array<double,9> box{{10,0,0,0,10,0,0,0,10}}, virial{};
    int step=17;
    p.cmd("setStep",&step);
    p.cmd("setPositions",positions.data());
    p.cmd("setMasses",masses.data());
    p.cmd("setBox",box.data());
    p.cmd("setForces",result.forces.data());
    p.cmd("setVirial",virial.data());
    p.cmd("calc");
  } catch(const PLMD::Plumed::Exception& error) {
    result.error=error.what();
  }
  return result;
}

bool localized(const Evaluation& result, const std::vector<std::string>& fields) {
  if(result.error.find("PLUMED action non-finite trace")==std::string::npos) {
    return false;
  }
  for(const auto& field : fields) {
    if(result.error.find(field)==std::string::npos) {
      return false;
    }
  }
  return true;
}

int main() {
  setenv("PLUMED_NONFINITE_ACTION_TRACE","d,c,r",1);
  std::ofstream output("output"), detail("details");
  bool passed=true;
  auto check=[&](bool success, const char* name, const Evaluation& result) {
    output << (success ? "PASS" : "FAIL") << " " << name << "\n";
    detail << name << "\n" << result.error << "\n";
    for(double force : result.forces) {
      detail << force << " ";
    }
    detail << "\n";
    passed=passed && success;
  };
  const Evaluation singular=evaluate({
    "d: DISTANCE ATOMS=1,2 NOPBC",
    "c: CUSTOM ARG=d FUNC=sqrt(x-1) PERIODIC=NO",
    "r: RESTRAINT ARG=c AT=1 KAPPA=1"}, {0,0,0,1,0,0});
  check(localized(singular, {"step=17","phase=forward","after_action=c ",
                             "target=c ","field=derivative"
                            }),
        "custom-derivative-localization",singular);

  const Evaluation zero=evaluate({
    "d: DISTANCE ATOMS=1,2 NOPBC",
    "r: RESTRAINT ARG=d AT=1 KAPPA=1"}, {0,0,0,0,0,0});
  check(localized(zero, {"phase=forward","after_action=d ",
                         "target=d ","field=derivative"
                        }),
        "colvar-derivative-localization",zero);

  // Finite CV, finite derivatives and finite scalar bias force can still
  // overflow when the chain rule is applied to atomic coordinates.
  const Evaluation overflow=evaluate({
    "c: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3 SELECT=1 "
    "KAPPA=1 REFERENCE=0.5 COEFFICIENTS=1e250 NOPBC",
    "r: RESTRAINT ARG=c AT=1 KAPPA=1e100"}, {0,0,0,2,0,0,1,0,0});
  check(localized(overflow, {"phase=backward","after_action=c ",
                             "target=posx ","field=force"
                            }),
        "atomic-force-overflow-localization",overflow);

  const Evaluation backward=evaluate({
    "d: DISTANCE ATOMS=1,2 NOPBC",
    "c: CUSTOM ARG=d FUNC=1e250*(x-1) PERIODIC=NO",
    "r: RESTRAINT ARG=c AT=1 KAPPA=1e100"}, {0,0,0,1,0,0});
  check(localized(backward, {"phase=backward","after_action=c ",
                             "target=d ","field=force"
                            }),
        "backward-force-localization",backward);

  const Evaluation boxOverflow=evaluate({
    "d: DISTANCE ATOMS=1,2 COMPONENTS NOPBC",
    "r: RESTRAINT ARG=d.x AT=1e250 KAPPA=0 SLOPE=1e100"},
  {0,0,0,1e250,0,0});
  check(localized(boxOverflow, {"phase=backward","after_action=d ",
                                "target=Box ","field=force"
                               }),
        "box-force-overflow-localization",boxOverflow);

  const Evaluation selected=evaluate({
    "d: DISTANCE ATOMS=1,2 NOPBC",
    "c: CUSTOM ARG=d FUNC=select(step(1.0-x),log(x+0.03),x-0.9704412) PERIODIC=NO",
    "r: RESTRAINT ARG=c AT=1 KAPPA=1"}, {0,0,0,1,0,0});
  bool finite=selected.error.empty();
  for(double force : selected.forces) {
    finite=finite && std::isfinite(force);
  }
  check(finite,"select-boundary-finite",selected);
  return passed ? 0 : 1;
}
