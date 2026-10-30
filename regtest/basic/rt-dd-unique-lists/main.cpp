#include "plumed/wrapper/Plumed.h"
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <string>
#include <vector>

// Serial test of the way DOMAIN_DECOMPOSITION shares positions and forces.
// The same calculation is done with three layouts of the arrays passed by the MD code:
//   0. all the atoms, in order
//   1. all the atoms, in reversed order (setAtomsGatindex)
//   2. only the atoms needed by PLUMED, as NAMD does (createFullList/getFullList)
// The input has actions that are only active at some steps (so that the atoms needed
// change from step to step), a neighbor list (so that the atoms requested by an action
// change during a step) and actions with and without forces.
// The test is run with PLUMED_FORCE_UNIQUE=auto, yes and no (see the config files of
// rt-dd-unique-lists*): the results must be identical in all cases.
// Finally, it checks that missing or too short arrays are reported as errors.

constexpr int natoms=100;
constexpr int nsteps=8;

typedef std::vector<std::array<double,3>> Forces;

static std::array<double,3> position(int atom,int step) {
  // atoms on a 5x5x4 grid in a box of side 2.5, with a small displacement that changes with time
  const int ix=atom%5, iy=(atom/5)%5, iz=atom/25;
  return {0.5*ix+0.05*std::sin(0.7*atom+0.3*step),
          0.5*iy+0.05*std::cos(1.3*atom+0.2*step),
          0.5*iz+0.05*std::sin(0.4*atom+0.5*step)};
}

static void readInput(PLMD::Plumed& p) {
  p.cmd("readInputLines",
        "d1: DISTANCE ATOMS=1,2\n"
        // the neighbor list is updated every 2 steps: at those steps all the atoms
        // are requested in prepare() and then only those in the list in calculate()
        "c1: COORDINATION GROUPA=3-10 GROUPB=11-30 R_0=0.3 NLIST NL_CUTOFF=0.8 NL_STRIDE=2\n"
        "RESTRAINT ARG=d1,c1 AT=0.4,2.0 KAPPA=10.0,1.0\n"
        // only active every 3 steps, no forces
        "t1: TORSION ATOMS=31,32,33,34\n"
        "PRINT ARG=t1 STRIDE=3 FILE=colvar_t1 FMT=%10.6f\n"
        // many atoms, only active every 4 steps: the atoms needed then exceed half of the atoms
        "c2: CENTER ATOMS=35-95\n"
        "d2: DISTANCE ATOMS=1,c2\n"
        "PRINT ARG=d2 STRIDE=4 FILE=colvar_d2 FMT=%10.6f\n"
        // only active every 5 steps, with forces
        "d3: DISTANCE ATOMS=96,97\n"
        "RESTRAINT ARG=d3 AT=0.4 KAPPA=5.0 STRIDE=5\n"
       );
}

// returns the forces on all the atoms at each step
static std::vector<Forces> run(int layout,std::vector<double>& biases) {
  PLMD::Plumed p;
  int n=natoms;
  double timestep=0.001;
  p.cmd("setNatoms",n);
  p.cmd("setMDEngine","test");
  p.cmd("setTimestep",timestep);
  p.cmd("setLogFile",("log"+std::to_string(layout)).c_str());
  p.cmd("init");
  readInput(p);
  std::vector<Forces> allforces;
  for(int step=0; step<nsteps; ++step) {
    p.cmd("setStep",step);
    // the global indexes of the atoms passed
    std::vector<int> list;
    if(layout==0) {
      for(int i=0; i<natoms; i++) {
        list.push_back(i);
      }
    } else if(layout==1) {
      for(int i=natoms-1; i>=0; i--) {
        list.push_back(i);
      }
    } else {
      p.cmd("prepareDependencies");
      int nlist=0;
      p.cmd("createFullList",&nlist);
      const int* plist=nullptr;
      p.cmd("getFullList",&plist);
      list.assign(plist,plist+nlist);
    }
    const int nlocal=list.size();
    std::vector<double> positions(3*nlocal), forces(3*nlocal,0.0), masses(nlocal,1.0), charges(nlocal,0.0);
    for(int i=0; i<nlocal; i++) {
      auto x=position(list[i],step);
      for(int k=0; k<3; k++) {
        positions[3*i+k]=x[k];
      }
    }
    double box[3][3]= {{2.5,0,0},{0,2.5,0},{0,0,2.0}}, virial[3][3]= {{0}};
    if(layout>0) {
      p.cmd("setAtomsNlocal",nlocal);
      p.cmd("setAtomsGatindex",list.data(), {static_cast<std::size_t>(nlocal)});
    }
    p.cmd("setPositions",positions.data(), {static_cast<std::size_t>(nlocal),3});
    p.cmd("setForces",forces.data(), {static_cast<std::size_t>(nlocal),3});
    p.cmd("setMasses",masses.data(), {static_cast<std::size_t>(nlocal)});
    p.cmd("setCharges",charges.data(), {static_cast<std::size_t>(nlocal)});
    p.cmd("setBox",&box[0][0], {3,3});
    p.cmd("setVirial",&virial[0][0], {3,3});
    if(layout==2) {
      p.cmd("shareData");
      p.cmd("performCalc");
    } else {
      p.cmd("calc");
    }
    double bias=0;
    p.cmd("getBias",&bias);
    biases.push_back(bias);
    Forces global(natoms, {0,0,0});
    for(int i=0; i<nlocal; i++) {
      for(int k=0; k<3; k++) {
        global[list[i]][k]=forces[3*i+k];
      }
    }
    allforces.push_back(global);
  }
  return allforces;
}

template<typename F>
static bool throws(F f, const std::string & expected) {
  try {
    f();
  } catch(const std::exception& e) {
    return std::string(e.what()).find(expected)!=std::string::npos;
  }
  return false;
}

int main(int argc,char**argv) {
  std::ofstream out("output");
  std::vector<std::vector<Forces>> results;
  std::vector<std::vector<double>> biases(3);
  for(int layout=0; layout<3; layout++) {
    results.push_back(run(layout,biases[layout]));
  }
  // the reference results are those of layout 0
  for(int step=0; step<nsteps; ++step) {
    char buffer[200];
    std::snprintf(buffer,200,"step %d bias %12.8f\n",step,biases[0][step]);
    out<<buffer;
    for(int i=0; i<natoms; i++) {
      const auto & f=results[0][step][i];
      if(f[0]!=0.0 || f[1]!=0.0 || f[2]!=0.0) {
        std::snprintf(buffer,200,"  atom %3d %12.8f %12.8f %12.8f\n",i+1,f[0],f[1],f[2]);
        out<<buffer;
      }
    }
  }
  // the other layouts must give the same results (up to the order of the operations)
  for(int layout=1; layout<3; layout++) {
    double maxdiff=0.0;
    for(int step=0; step<nsteps; ++step) {
      maxdiff=std::max(maxdiff,std::fabs(biases[layout][step]-biases[0][step]));
      for(int i=0; i<natoms; i++)
        for(int k=0; k<3; k++) {
          maxdiff=std::max(maxdiff,std::fabs(results[layout][step][i][k]-results[0][step][i][k]));
        }
    }
    out<<"layout "<<layout<<" agrees with layout 0: "<<(maxdiff<1e-12?"yes":"NO")<<std::endl;
  }
  // with an argument, the checks of the errors are skipped (to compare with older versions)
  if(argc>1) {
    return 0;
  }

  // errors: forces applied on the atoms, but no array for the forces
  out<<"missing forces array detected: "<<(throws([]() {
    PLMD::Plumed p;
    int n=natoms;
    p.cmd("setNatoms",n);
    p.cmd("setLogFile","/dev/null");
    p.cmd("init");
    p.cmd("readInputLines","d: DISTANCE ATOMS=1,2\nRESTRAINT ARG=d AT=0 KAPPA=1\n");
    std::vector<double> positions(3*natoms,0.0), masses(natoms,1.0);
    positions[3]=1.0;
    double box[3][3]= {{2.5,0,0},{0,2.5,0},{0,0,2.0}}, virial[3][3]= {{0}};
    int step=0;
    p.cmd("setStep",step);
    p.cmd("setPositions",positions.data(), {static_cast<std::size_t>(natoms),3});
    p.cmd("setMasses",masses.data(), {static_cast<std::size_t>(natoms)});
    p.cmd("setBox",&box[0][0], {3,3});
    p.cmd("setVirial",&virial[0][0], {3,3});
    p.cmd("calc");
  },"did not pass the array")?"yes":"NO")<<std::endl;

  // errors: a bias on a value passed with PUT, but no array for its force
  out<<"missing value force detected: "<<(throws([]() {
    PLMD::Plumed p;
    int n=0;
    p.cmd("setNatoms",n);
    p.cmd("setLogFile","/dev/null");
    p.cmd("init");
    p.cmd("readInputLines","ext: PUT UNIT=number PERIODIC=NO MUTABLE\nRESTRAINT ARG=ext AT=0 KAPPA=1\n");
    double val=1.0;
    int step=0;
    p.cmd("setStep",step);
    p.cmd("setValue ext",&val);
    p.cmd("calc");
  },"did not pass the array")?"yes":"NO")<<std::endl;

  // errors: the positions of shuffled atoms are passed in an array that is too short
  out<<"short positions array detected: "<<(throws([]() {
    PLMD::Plumed p;
    int n=natoms;
    p.cmd("setNatoms",n);
    p.cmd("setLogFile","/dev/null");
    p.cmd("init");
    p.cmd("readInputLines","d: DISTANCE ATOMS=1,2\nPRINT ARG=d FILE=/dev/null\n");
    std::vector<int> list;
    for(int i=natoms-1; i>=0; i--) {
      list.push_back(i);
    }
    std::vector<double> positions(3*natoms,0.0), forces(3*natoms,0.0), masses(natoms,1.0);
    double box[3][3]= {{2.5,0,0},{0,2.5,0},{0,0,2.0}}, virial[3][3]= {{0}};
    for(int step=0; step<2; step++) {
      p.cmd("setStep",step);
      p.cmd("setAtomsNlocal",n);
      p.cmd("setAtomsGatindex",list.data(), {static_cast<std::size_t>(natoms)});
      // one atom less than the local atoms
      p.cmd("setPositions",positions.data(), {static_cast<std::size_t>(natoms-1),3});
      p.cmd("setForces",forces.data(), {static_cast<std::size_t>(natoms),3});
      p.cmd("setMasses",masses.data(), {static_cast<std::size_t>(natoms)});
      p.cmd("setBox",&box[0][0], {3,3});
      p.cmd("setVirial",&virial[0][0], {3,3});
      p.cmd("calc");
    }
  },"elements on axis 0")?"yes":"NO")<<std::endl;
  return 0;
}
