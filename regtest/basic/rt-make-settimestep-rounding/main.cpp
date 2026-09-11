#include <vector>
#include <fstream>
#include "plumed/wrapper/Plumed.h"
#include <iostream>

using namespace PLMD;

template<typename T>
void run(){

  auto natoms=10;
  std::vector<T> masses;
  std::vector<std::array<T,3>> positions;
  std::vector<std::array<T,3>> forces;
  T box[3][3];
  T virial[3][3];
  
  Plumed p;

  p.cmd("setRealPrecision",int(sizeof(T)));
  p.cmd("setNatoms",natoms);
  p.cmd("setTimestep",(T)0.002);
  p.cmd("init");
  
  p.cmd("readInputLines",
    "c: CONSTANT VALUE=1.0 \n"
    "PRINT ARG=c FILE=COLVAR RESTART=YES\n"
  );

  // dummy settings, not really used
  positions.resize(natoms);
  masses.resize(natoms);
  forces.resize(natoms);

  for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++) positions[i][j]=0.0;
  for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++) forces[i][j]=0.0;
  for(unsigned i=0;i<natoms;i++) masses[i]=1.0;

  p.cmd("setStep",10000000);
  p.cmd("setMasses",&masses[0]);
  p.cmd("setPositions",&positions[0][0]);
  p.cmd("setForces",&forces[0][0]);

  p.cmd("setBox",&box[0][0]);
  p.cmd("setVirial",&virial[0][0]);
  p.cmd("calc");

}

int main(){
  run<double>();
  run<float>();
}
