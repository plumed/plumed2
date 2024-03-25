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
    throw std::runtime_error("this should have failed");
  } catch(const Plumed::ExceptionTypeError & err) {
    // expected
  }
  p.cmd("setStep",0);
  //p.cmd("setPositions",&positions[0][0],{natoms,3}); // not c++98
  std::size_t shape[3] = {natoms,3,0};
  p.cmd("setPositions",&positions[0][0],shape);

  p.cmd("setVirial",&virial[0][0]);
  p.cmd("setBox",&box[0][0]);
  p.cmd("setForces",&forces[0][0]);
  p.cmd("setMasses",&masses[0]);

  p.cmd("calc");
  
  return 0;
}
