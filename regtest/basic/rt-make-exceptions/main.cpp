#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/Exception.h"
#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace PLMD;

void test_line(std::ostream & ofs,Plumed & p,const std::string & arg){
  std::string cmd="readInputLine";
  ofs<<cmd<<" "<<arg<<std::endl;
  try{
    p.cmd(cmd.c_str(),arg.c_str());
    ofs<<"+++ !!!! uncatched !!!!"<<std::endl;
  } catch(Exception&e) {
    ofs<<"+++ catched"<<std::endl;
  }
}

void test_this(std::ostream & ofs,Plumed & p,const std::string & cmd,const void*arg){
  ofs<<cmd;
  if(!arg) ofs<<" NULL";
  ofs<<std::endl;
  try{
    p.cmd(cmd.c_str(),arg);
    ofs<<"+++ !!!! uncatched !!!!"<<std::endl;
  } catch(Exception&e) {
    ofs<<"+++ catched"<<std::endl;
  }
}

int main(){

  std::ofstream ofs("output");

  {
// test a mistake in timer
    Stopwatch sw;
    ofs<<"pause"<<std::endl;
    try{
      sw.pause();
      ofs<<"+++ !!!! uncatched !!!!"<<std::endl;
    } catch(Exception& e) {
      ofs<<"+++ catched"<<std::endl;
    }
  }

  Plumed plumed;

// first try to wrongly set real precision
  unsigned i=127;
  test_this(ofs,plumed,"setRealPrecision",&i);

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
  test_line(ofs,plumed,"d: DISTANCE ATOMS=1,2,3");
  test_line(ofs,plumed,"d:DISTANCE ATOMS=1,2");
  test_line(ofs,plumed,"d: DIST ANCE ATOMS=1,2");
  test_line(ofs,plumed,"d: DISTANCE ATOMS=1,2 COMPONENTS SCALED_COMPONENTS");
  test_line(ofs,plumed,"d: GYRATION ATOMS=");
  test_line(ofs,plumed,"d: GYRATION ATOMS=1-4 TYPE=WHAT");
  test_line(ofs,plumed,"d: POSITION ATOM=1,2");
  test_line(ofs,plumed,"d: PUCKERING ATOMS=1-4");
  test_line(ofs,plumed,"d: ANGLE ATOMS=1,2,3,4,5");
  test_line(ofs,plumed,"d: COORDINATION GROUPA=1 GROUPB=2 R_0=0.5 NN=1.5");

// these should not fail
  plumed.cmd("readInputLine","d: DISTANCE ATOMS=1,2");
  plumed.cmd("readInputLine","d1: DISTANCE ATOMS={1 2}"); // check if braces are parsed correctly
  plumed.cmd("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1");

  test_line(ofs,plumed,"METAD ARG=d PACE=1 SIGMA=1 HEIGHT=0 FILE=H1 RESTART=WHAT");
  test_line(ofs,plumed,"METAD ARG=d PACE=1 SIGMA=1 TAU=5");
  test_line(ofs,plumed,"COMBINE ARG=d,d1 COEFFICIENTS=3");
  test_line(ofs,plumed,"COMBINE ARG=d,d1 COEFFICIENTS=3,3 PARAMETERS=1");
  test_line(ofs,plumed,"COMBINE ARG=d,d1 COEFFICIENTS=3,3 PARAMETERS=1,2 POWERS=4");

// these should not fail
  plumed.cmd("readInputLine","m1: METAD ARG=d PACE=1 SIGMA=5 HEIGHT=1 FILE=H1 FMT=%9.5f");
  plumed.cmd("readInputLine","m2: METAD ARG=d PACE=2 SIGMA=5 HEIGHT=1 FILE=H2 FMT=%9.5f");
  plumed.cmd("readInputLine","PRINT ARG=d,d1,m1.bias FILE=COLVAR FMT=%9.5f");

  test_this(ofs,plumed,"something random here",NULL);
  for(int step=0;step<3;step++){

// this should fail
    test_this(ofs,plumed,"setStep",NULL);

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
    test_this(ofs,plumed,"setMasses",&masses[0]);
  }

  return 0;

}
