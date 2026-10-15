#include "plumed/wrapper/Plumed.h"

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct TestCase {
  std::string name;
  std::string input;
  std::vector<double> positions;
  std::array<double,9> box;
  std::vector<std::string> expected;
  bool shouldFail;
};

std::string evaluate(const TestCase& test) {
  try {
    PLMD::Plumed plumed;
    int precision=8;
    int natoms=static_cast<int>(test.positions.size()/3);
    double timestep=0.001;
    plumed.cmd("setRealPrecision",&precision);
    plumed.cmd("setNatoms",&natoms);
    plumed.cmd("setTimestep",&timestep);
    plumed.cmd("init");
    plumed.cmd("readInputLine",test.input.c_str());
    plumed.cmd("readInputLine","PRINT ARG=v FILE=/dev/null");

    std::vector<double> positions=test.positions;
    std::vector<double> forces(positions.size(),0.0);
    std::vector<double> masses(static_cast<unsigned>(natoms),1.0);
    std::array<double,9> box=test.box;
    std::array<double,9> virial{{0,0,0,0,0,0,0,0,0}};
    int step=7;
    plumed.cmd("setStep",&step);
    plumed.cmd("setPositions",positions.data());
    plumed.cmd("setMasses",masses.data());
    plumed.cmd("setBox",box.data());
    plumed.cmd("setForces",forces.data());
    plumed.cmd("setVirial",virial.data());
    plumed.cmd("calc");
  } catch(const PLMD::Plumed::Exception& error) {
    return error.what();
  }
  return std::string();
}

bool checkCase(std::ofstream& output, const TestCase& test) {
  const std::string message=evaluate(test);
  if(test.shouldFail && message.empty()) {
    output << "FAIL " << test.name << " did not fail\n";
    return false;
  }
  if(!test.shouldFail && !message.empty()) {
    output << "FAIL " << test.name << " unexpectedly failed\n";
    return false;
  }
  for(const std::string& expected : test.expected) {
    if(message.find(expected)==std::string::npos) {
      output << "FAIL " << test.name << " missing " << expected << "\n";
      return false;
    }
  }
  output << "PASS " << test.name << "\n";
  return true;
}

std::array<double,9> diagonalBox(const double x, const double y,
                                 const double z) {
  return std::array<double,9> {{x,0,0,0,y,0,0,0,z}};
}

}

int main() {
  const double epsilon=std::numeric_limits<double>::epsilon();
  const double nan=std::numeric_limits<double>::quiet_NaN();
  const double infinity=std::numeric_limits<double>::infinity();
  const std::string normalInput=
    "v: VORONOI_COORDINATION CENTERS=1 ASSIGNED=2 KAPPA=5 "
    "REFERENCE=1 POWER=2 NOPBC";
  const std::string pbcInput=
    "v: VORONOI_COORDINATION CENTERS=1 ASSIGNED=2 KAPPA=5 "
    "REFERENCE=1 POWER=2";
  const std::array<double,9> box=diagonalBox(3.0,3.0,3.0);
  std::vector<TestCase> tests;
  tests.push_back({"normal",normalInput,{0,0,0,1,0,0},box,{},false});
  tests.push_back({"exact-overlap",normalInput,{0,0,0,0,0,0},box,
    {
      "a CENTER-ASSIGNED distance is zero or non-finite",
      "category=distance-zero-or-near-zero","step=7","pair_index=0(zero-based)",
      "center_atom=1 assigned_atom=2","coordinate_units=PLUMED-internal"
    },true});
  tests.push_back({"pbc-image-overlap",pbcInput,{0,0,0,3,0,0},box,
    {
      "category=distance-zero-or-near-zero","minimum_image_displacement=(0,0,0)",
      "pbc=on","box_determinant=27"
    },true});
  tests.push_back({"below-epsilon",normalInput,{0,0,0,0.5*epsilon,0,0},box,
    {"category=distance-zero-or-near-zero"},true});
  tests.push_back({"above-epsilon",normalInput,{0,0,0,2.0*epsilon,0,0},box,{},false});
  tests.push_back({"center-nan",normalInput,{nan,0,0,1,0,0},box,
    {"category=center-coordinate-non-finite","center_position=(nan,0,0)"},true});
  tests.push_back({"assigned-inf",normalInput,{0,0,0,infinity,0,0},box,
    {"category=assigned-coordinate-non-finite","assigned_position=(inf,0,0)"},true});
  tests.push_back({"displacement-overflow",normalInput,
    {-1.7e308,0,0,1.7e308,0,0},box,
    {"category=minimum-image-displacement-non-finite"},true});
  tests.push_back({"length-overflow",normalInput,{0,0,0,1.0e200,0,0},box,
    {"category=distance-non-finite","length=inf"},true});
  tests.push_back({"score-overflow",
                   "v: VORONOI_COORDINATION CENTERS=1 ASSIGNED=2 KAPPA=1e308 "
  "REFERENCE=1 POWER=2 NOPBC",{0,0,0,2,0,0},box,
  {
    "KAPPA times a CENTER-ASSIGNED distance is too large",
    "category=kappa-times-distance-non-finite"
  },true});
  tests.push_back({"degenerate-box",pbcInput,{0,0,0,1,0,0},
                   diagonalBox(0.0,3.0,3.0),
  {"category=pbc-box-degenerate","pbc=on","box_determinant=0"},true});

  std::vector<double> manyPositions(32*3,0.0);
  for(unsigned i=0; i<16; ++i) {
    manyPositions[3*i]=static_cast<double>(i);
    manyPositions[3*(16+i)]=100.0+static_cast<double>(i);
  }
  manyPositions[3*31]=15.0;
  tests.push_back({"openmp-deterministic-first-pair",
                   "v: VORONOI_COORDINATION CENTERS=1-16 ASSIGNED=17-32 KAPPA=5 "
                   "REFERENCE=1 POWER=2 NOPBC",manyPositions,diagonalBox(300.0,300.0,300.0),
  {
    "category=distance-zero-or-near-zero","pair_index=255(zero-based)",
    "center_atom=16 assigned_atom=32","rank=0"
  },true});

  std::ofstream output("output");
  bool passed=true;
  for(const TestCase& test : tests) {
    passed=checkCase(output,test) && passed;
  }
  return passed ? 0 : 1;
}
