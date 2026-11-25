#include "type_traits"
// custom structure used to store 3d vectors
// this is similar to PLMD::Vector and to OpenMM::Vec3
class vec3 {
  double data[3];
public:
  double & operator[](unsigned i) { return data[i]; }
  const double & operator[](unsigned i) const { return data[i]; }
};

// users need to add this so as to make plumed aware of their structure
namespace PLMD {
namespace wrapper {
// forward declaration
template<typename T> struct is_custom_array;

template<> struct is_custom_array<vec3> : std::true_type {
  // this is just saying that the elements of this array are double
  // plumed is going to compute the number of elements as sizeof(vec3)/sizeof(double)
  using value_type = double;
};
}
}

// In strict mode:
// - full shapes should be passed to plumed
// - when passing pointers with no shape information, shape is assumed to be {1}
#define __PLUMED_WRAPPER_CXX_DETECT_SHAPES_STRICT 1
#include "plumed/wrapper/Plumed.h"
#include <fstream>

using namespace PLMD;



int main(){

  unsigned natoms=10;

  std::vector<double> masses;
  std::vector<vec3> positions;
  std::vector<vec3> forces;
  vec3 box[3];
  vec3 virial[3];
  double bias;

  Plumed p;

  p.cmd("setNatoms",natoms);
  p.cmd("init");

  p.cmd("readInputLines",
    "d: DISTANCE ATOMS=1,10\n"
    "RESTRAINT ARG=d AT=0 KAPPA=10\n"
    "c: CELL\n"
    "PRINT ARG=d,c.* FILE=COLVAR\n"
  );

  positions.resize(natoms);
  masses.resize(natoms);
  forces.resize(natoms);

  for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++) positions[i][j]=i*10+j;
  for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++) forces[i][j]=0.0;
  for(unsigned i=0;i<natoms;i++) masses[i]=i;

  box[0][0]=10.0; box[0][1]=0.0; box[0][2]=0.0; box[1][0]=1.0; box[1][1]=10.0; box[1][2]=0.0; box[2][0]=2.0; box[2][1]=1.0; box[2][2]=10.0;
  for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) virial[i][j]=0.0;

  p.cmd("setStep",0);

  //// this is the best option:
  p.cmd("setMasses",masses);

  //// this would be fine, but it's redundant
  // p.cmd("setMasses",masses.data(),{natoms});

  //// this would trigger a compilation error, because it is not a full shape
  //// it is instead allowed when not using the strict mode
  // p.cmd("setMasses",masses.data(),natoms);

  //// this would trigger a runtime error later (in cmd("calc")), because it is equivalent to p.cmd("setMasses",{1})
  //// it is instead allowed when not using the strict mode
  // p.cmd("setMasses",masses.data());

  //// this is the best option
  p.cmd("setPositions",positions);

  //// this would be fine, but it's redundant
  // p.cmd("setPositions",&positions[0][0],{natoms,3});

  //// this would trigger a compilation error, because it is not a full shape
  //// it is instead allowed when not using the strict mode
  // p.cmd("setPositions",&positions[0][0],3*natoms);

  //// this would trigger a runtime error later (in cmd("calc")), because it is equivalent to p.cmd("setPositions",{1})
  //// it is instead allowed when not using the strict mode
  // p.cmd("setPositions",&positions[0][0]);

  p.cmd("setForces",forces);

  //// this was not possible before because box would have been converted to void*
  //// now, this is the recommended options and notifies plumed that this ia a 3x3 tensor
  p.cmd("setBox",box);

  p.cmd("setVirial",virial);
  p.cmd("calc");

  //// this is ok, and is equivalent to p.cmd("getBias",&bias,{1})
  p.cmd("getBias",&bias);


  // write the results on a file
  std::ofstream os("ooo");
  for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++) os<<" "<<forces[i][j];
  os<<"\n";
  for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) os<<" "<<virial[i][j];
  os<<"\n";
  os<<bias<<"\n";

  return 0;
  
}
