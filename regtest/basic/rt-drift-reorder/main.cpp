#include "plumed/wrapper/Plumed.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Serial test of EFFECTIVE_ENERGY_DRIFT when the MD code changes the order of the
// atoms from one step to the next (as LAMMPS does when it sorts the atoms).
// The same trajectory is passed with the atoms in order (layout 0) and with
// a different order at every step (layout 1): the effective energy drift must be the same.

constexpr int natoms=20;
constexpr int nsteps=12;

static double coordinate(int atom,int k,int step) {
  // a slow motion, crossing the boundaries of the box for some atoms
  return 0.1*(atom%4+k) + 0.9*k + 0.02*step*(1+(atom+k)%3) + 0.05*std::sin(0.7*atom+0.4*step+k);
}

static std::vector<double> run(int layout) {
  std::string eff="eff"+std::to_string(layout);
  {
  PLMD::Plumed p;
  int n=natoms;
  double timestep=0.002;
  p.cmd("setNatoms",n);
  p.cmd("setMDEngine","test");
  p.cmd("setTimestep",timestep);
  p.cmd("setLogFile",("log"+std::to_string(layout)).c_str());
  p.cmd("init");
  p.cmd("readInputLines",
        ("d1: DISTANCE ATOMS=1,10\n"
         "RESTRAINT ARG=d1 AT=0.5 KAPPA=100\n"
         "c: CENTER ATOMS=1-20 NOPBC\n"
         "d2: DISTANCE ATOMS=c,5\n"
         "RESTRAINT ARG=d2 AT=0.2 KAPPA=50\n"
         "EFFECTIVE_ENERGY_DRIFT PRINT_STRIDE=1 FILE="+eff+" FMT=%14.10f\n").c_str());
  for(int step=0; step<nsteps; ++step) {
    p.cmd("setStep",step);
    std::vector<int> list(natoms);
    for(int i=0; i<natoms; i++) {
      // layout 1: a rotation that changes at every step, reversed at odd steps
      int j=(i+3*step)%natoms;
      list[i]=(layout==0 ? i : (step%2 ? natoms-1-j : j));
    }
    std::vector<double> positions(3*natoms), forces(3*natoms,0.0), masses(natoms,1.0), charges(natoms,0.0);
    for(int i=0; i<natoms; i++)
      for(int k=0; k<3; k++) {
        positions[3*i+k]=coordinate(list[i],k,step);
      }
    double box[3][3]= {{2.0,0,0},{0,2.2,0},{0,0,2.4}}, virial[3][3]= {{0}};
    if(layout>0) {
      p.cmd("setAtomsNlocal",n);
      p.cmd("setAtomsGatindex",list.data(), {static_cast<std::size_t>(natoms)});
    }
    p.cmd("setPositions",positions.data(), {static_cast<std::size_t>(natoms),3});
    p.cmd("setForces",forces.data(), {static_cast<std::size_t>(natoms),3});
    p.cmd("setMasses",masses.data(), {static_cast<std::size_t>(natoms)});
    p.cmd("setCharges",charges.data(), {static_cast<std::size_t>(natoms)});
    p.cmd("setBox",&box[0][0], {3,3});
    p.cmd("setVirial",&virial[0][0], {3,3});
    p.cmd("calc");
  }
  } // the file is closed when the Plumed object is destroyed
  // read back the effective energy drift
  std::vector<double> result;
  std::ifstream in(eff);
  std::string line;
  while(std::getline(in,line)) {
    if(line.empty() || line[0]=='#') {
      continue;
    }
    std::istringstream is(line);
    double t,e;
    is>>t>>e;
    result.push_back(e);
  }
  return result;
}

int main() {
  std::ofstream out("output");
  auto ref=run(0);
  auto shuffled=run(1);
  char buffer[200];
  for(unsigned i=0; i<ref.size(); i++) {
    std::snprintf(buffer,200,"step %u eed %14.10f\n",i,ref[i]);
    out<<buffer;
  }
  double maxdiff=(ref.size()==shuffled.size() && ref.size()==nsteps) ? 0.0 : 1.0;
  for(unsigned i=0; i<ref.size() && i<shuffled.size(); i++) {
    maxdiff=std::max(maxdiff,std::fabs(ref[i]-shuffled[i]));
  }
  out<<"shuffled atoms give the same drift: "<<(maxdiff<1e-9?"yes":"NO")<<std::endl;
  return 0;
}
