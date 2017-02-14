#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/Exception.h"
#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace PLMD;

int main(){

  std::ofstream ofs("output");

  {
// test a mistake in timer
    Stopwatch sw;
    try{ sw.pause();
    } catch(Exception& e) { ofs<<"E pause\n"; }
  }

  Plumed plumed;

// first try to wrongly set real precision
  unsigned i=127;
  try{ plumed.cmd("setRealPrecision",&i);
  } catch(Exception&e){ ofs<<"E setRealPrecision"<<std::endl;}

  int natoms=10;

  std::vector<double> positions(3*natoms,0.0);
  for(unsigned i=0;i<3*natoms;i++) positions[i]=i;
  std::vector<double> masses(natoms,1.0);
  std::vector<double> forces(3*natoms,0.0);
  std::vector<double> box(9,0.0);
  std::vector<double> virial(9,0.0);

  plumed.cmd("setNatoms",&natoms);
  plumed.cmd("setLogFile","test.log");
  plumed.cmd("init");

// I try many mistaken lines.
// Each of them will raise an exception
// Notice that name "d" will not be reserved and it will be possible
// to use it later
  try{ plumed.cmd("readInputLine","d: DISTANCE ATOMS=1,2,3");
  } catch(Exception& e) {   ofs<<"E d:DISTANCE ATOMS=1,2,3"<<std::endl; }
  try{ plumed.cmd("readInputLine","d:DISTANCE ATOMS=1,2");
  } catch(Exception& e) {   ofs<<"E d:DISTANCE ATOMS=1,2"<<std::endl; }
  try{ plumed.cmd("readInputLine","d: DIST ANCE ATOMS=1,2");
  } catch(Exception& e) {   ofs<<"E d: DIST ANCE ATOMS=1,2"<<std::endl; }
  try{ plumed.cmd("readInputLine","d: ANGLE ATOMS=1,2,3,4,5");
  } catch(Exception& e) {   ofs<<"E d: ANGLE ATOMS=1,2,3,4,5 "<<std::endl; }
  try{ plumed.cmd("readInputLine","d: COORDINATION GROUPA=1 GROUPB=2 R_0=0.5 NN=1.5");
  } catch(Exception& e) {   ofs<<"E d: COORDINATION GROUPA=1 GROUPB=2 R_0=0.5 NN=1.5"<<std::endl; }

// these should not fail
  plumed.cmd("readInputLine","d: DISTANCE ATOMS=1,2");
  plumed.cmd("readInputLine","d1: DISTANCE ATOMS={1 2}"); // check if braces are parsed correctly
  plumed.cmd("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1");

// Check stupid option to RESTART
  try{ plumed.cmd("readInputLine","METAD ARG=d PACE=1 SIGMA=1 HEIGHT=0 FILE=H1 RESTART=WHAT");
  } catch(Exception& e) {   ofs<<"E METAD ARG=d PACE=1 SIGMA=1 HEIGHT=0 FILE=H1 RESTART=WHAT"<<std::endl; }

// these should not fail
  plumed.cmd("readInputLine","m1: METAD ARG=d PACE=1 SIGMA=5 HEIGHT=1 FILE=H1 FMT=%9.5f");
  plumed.cmd("readInputLine","m2: METAD ARG=d PACE=2 SIGMA=5 HEIGHT=1 FILE=H2 FMT=%9.5f");
  plumed.cmd("readInputLine","PRINT ARG=d,d1,m1.bias FILE=COLVAR FMT=%9.5f");

  try{ plumed.cmd("something random here",NULL);
  } catch(Exception& e) { ofs<<"E random cmd"<<std::endl;}
  for(int step=0;step<3;step++){

// this should fail
    try{ plumed.cmd("setStep",NULL);
    } catch(Exception& e) { ofs<<"E cmd setStep NULL"<<std::endl;}

    plumed.cmd("setStep",&step);
    plumed.cmd("setPositions",&positions[0]);
    plumed.cmd("setBox",&box[0]);
    plumed.cmd("setForces",&forces[0]);
    plumed.cmd("setVirial",&virial[0]);
    plumed.cmd("setMasses",&masses[0]);
// set positions after having passed the pointer. They should be accessed here (at "calc").
    for(unsigned i=0;i<3*natoms;i++) positions[i]=i*step;
    plumed.cmd("calc");

// this should fail
    try{ plumed.cmd("setMasses",&masses[0]);
    } catch(Exception& e) { ofs<<"E setMasses called in wrong place"<<std::endl;}
  }

  return 0;

}
