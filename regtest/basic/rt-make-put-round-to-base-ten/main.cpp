#include <vector>
#include <array>
#include "plumed/wrapper/Plumed.h"

using namespace PLMD;

// Pass two scalars from a single precision MD code: one that the user writes in
// base ten, such as the timestep, and one that the MD code computes, such as the
// potential energy.  Only the former should be rounded, and a negative value
// should not become nan.
template<typename T>
void run() {

  unsigned natoms=10;
  std::vector<T> masses(natoms,1.0);
  std::vector<std::array<T,3>> positions(natoms);
  std::vector<std::array<T,3>> forces(natoms);
  T box[3][3]= {{0.0}};
  T virial[3][3]= {{0.0}};

  T dt=0.002;
  T energy=-1234.5678;

  Plumed p;
  p.cmd("setRealPrecision",int(sizeof(T)));
  p.cmd("setNatoms",natoms);
  p.cmd("init");

  p.cmd("readInputLines",
        "dt_plain: PUT UNIT=number PERIODIC=NO CONSTANT SHAPE=0\n"
        "dt_round: PUT UNIT=number PERIODIC=NO CONSTANT ROUND_TO_BASE_TEN SHAPE=0\n"
        "en_plain: PUT UNIT=number PERIODIC=NO CONSTANT SHAPE=0\n"
        "en_round: PUT UNIT=number PERIODIC=NO CONSTANT ROUND_TO_BASE_TEN SHAPE=0\n"
        "PRINT ARG=dt_plain,dt_round,en_plain,en_round FILE=COLVAR FMT=%.12f RESTART=YES\n"
       );

  p.cmd("setValue dt_plain",&dt);
  p.cmd("setValue dt_round",&dt);
  p.cmd("setValue en_plain",&energy);
  p.cmd("setValue en_round",&energy);

  for(unsigned i=0; i<natoms; i++) {
    for(unsigned j=0; j<3; j++) {
      positions[i][j]=0.0;
      forces[i][j]=0.0;
    }
  }

  p.cmd("setStep",1);
  p.cmd("setMasses",&masses[0]);
  p.cmd("setPositions",&positions[0][0]);
  p.cmd("setForces",&forces[0][0]);
  p.cmd("setBox",&box[0][0]);
  p.cmd("setVirial",&virial[0][0]);
  p.cmd("calc");
}

int main() {
  run<double>();
  run<float>();
}
