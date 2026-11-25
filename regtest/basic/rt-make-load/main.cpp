#include "plumed/wrapper/Plumed.h"
#include <vector>
#include <thread>

using namespace PLMD;

int natoms=4;

void init(Plumed &p,std::string logfile) {
  p.cmd("setNatoms",natoms);
  p.cmd("setLogFile",logfile.c_str());
  p.cmd("init");
}

void readline(Plumed &p,const std::string & line) {
  p.cmd("readInputLine",line.c_str());
}

void run(Plumed &p) {
  std::vector<std::array<double,3>> positions(natoms);
  std::vector<std::array<double,3>> forces(natoms);
  std::vector<double> masses(natoms);
  for(auto i=0; i<natoms; i++) masses[i]=i+1;
  for(auto i=0; i<natoms; i++) for(auto j=0;j<3;j++) positions[i][j]=10*i+j;
  double cell[3][3];
  double virial[3][3];
  for(auto i=0; i<3; i++) for(auto j=0;j<3;j++) cell[i][j]=0.0;
  for(auto i=0; i<3; i++) for(auto j=0;j<3;j++) virial[i][j]=0.0;
  cell[0][0]=100.0;
  cell[1][1]=100.0;
  cell[2][2]=100.0;
  p.cmd("setStep",0);
  p.cmd("setBox",&cell[0][0]);
  p.cmd("setVirial",&virial[0][0]);
  p.cmd("setMasses",masses.data());
  p.cmd("setPositions",&positions[0][0]);
  p.cmd("setForces",&forces[0][0]);
  p.cmd("calc");
}

// i is thread number (starting from 0)
void fun(Plumed p,unsigned i) {
  init(p,std::string("log_threads") + std::to_string(i));
  
  readline(p,"d1: DISTANCE ATOMS=1,2");
  // threads 0,2,4 etc are loading Distance10
  // threads 1,3,5 etc are loading Distance 20
  readline(p,"LOAD FILE=./Distance" + std::to_string((i%2+1)*10) + "." __PLUMED_SOEXT);
  readline(p,"d_loaded: DISTANCE ATOMS=1,2");
  // output is written on a separate file for each thread
  readline(p,"PRINT FILE=test_threads" + std::to_string(i) + " ARG=d1,d_loaded");
  run(p);
}

int main() {
  // threads
  {
    unsigned nthreads=16;
    std::vector<std::thread> threads;
    // plumed objects are declared outside
    std::vector<Plumed> pl(nthreads);
    for(unsigned j=0;j<nthreads;j++) threads.emplace_back(fun,pl[j],j);
    // threads are joined before destruction
    for(unsigned j=0;j<nthreads;j++) threads[j].join();

    // this is done to make sure that all plumed objects survive till the end
    // and make the number of unique dylib loaded reproducible
  }

  Plumed p1,p2;

  init(p1,"log_sequential1");
  init(p2,"log_sequential2");

  // sequential, controlling the loading order
  readline(p1,"d1: DISTANCE ATOMS=1,2");
  readline(p1,"LOAD FILE=./Distance10." __PLUMED_SOEXT);
  readline(p1,"d10: DISTANCE ATOMS=1,2");
  readline(p1,"PRINT FILE=test1 ARG=d1,d10");
  readline(p2,"d1: DISTANCE ATOMS=1,2");
  readline(p2,"LOAD FILE=./Distance20." __PLUMED_SOEXT);
  readline(p2,"d20: DISTANCE ATOMS=1,2");
  readline(p2,"PRINT FILE=test2 ARG=d1,d20");

  run(p1);
  run(p2);


  return 0;
}

