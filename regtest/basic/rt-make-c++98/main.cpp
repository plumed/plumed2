/*
This test is useful to check what happens if Plumed.h is included in a program
not using C++11. This is a very basic test, and it is here mostly to verify
that the various #ifdefs in Plumed.h do not result in something syntactically wrong
when using C++ pre 11.
*/
#include "plumed/wrapper/Plumed.h"
#include <vector>

using namespace PLMD;

class vec3 {
  double pos[3];
public:
  double & operator[](unsigned i){ return pos[i];}
  const double & operator[](unsigned i) const { return pos[i];}
};

int main() {
  PLMD::Plumed p;
  unsigned natoms=10;
  std::vector<vec3> positions(natoms);
  std::vector<vec3> forces(natoms);
  std::vector<double> masses(natoms);
  for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++) positions[i][j]=i+j;
  double virial[3][3];
  double box[3][3];
  for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) box[i][j]=0.0;

  box[0][0]=box[1][1]=box[2][2]=10.0;

  p.cmd("setNatoms",natoms);
  p.cmd("init");
  p.cmd("readInputLine","d: GYRATION ATOMS=1-9");
  p.cmd("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1");
  p.cmd("readInputLine","PRINT FILE=COLVAR ARG=d RESTART=YES");
  try {
    p.cmd("setStep",0.0);
    throw std::runtime_error("this should have failed (argument is float instead of int)");
  } catch(const Plumed::ExceptionTypeError & err) {
    // expected
  }
  p.cmd("setStep",0);

  // this is not allowed in c++98:
  //p.cmd("setPositions",&positions[0][0],{natoms,3});
  // we should rather use an explicitly null terminated array:
  std::size_t shape[3] = {natoms,3,0};
  p.cmd("setPositions",&positions[0][0],shape);

  p.cmd("setVirial",&virial[0][0]);
  p.cmd("setBox",&box[0][0]);
  p.cmd("setForces",&forces[0][0]);
  p.cmd("setMasses",&masses[0]);

  p.cmd("calc");
  
  try {
    p.cmd("setStep",0);
    p.cmd("setPositions",&positions[0][0]);
    p.cmd("setVirial",&virial[0][0]);
    p.cmd("setBox",box);
    p.cmd("setForces",&forces[0][0]);
    p.cmd("setMasses",&masses[0]);
    p.cmd("calc");
    throw std::runtime_error("this should have failed (argument is void* instead of double*)");
  } catch(const Plumed::ExceptionTypeError & err) {
    // expected, since in C++98 we cannot enable shape detection, and box is seen as a void*.
    // the same thing was actually happening with PLUMED 2.9, even when using a C++11 compiler,
    // but should not happen in PLUMED 2.10, when using a C++11 compiler
  }

  return 0;
}
