#include "plumed/wrapper/Plumed.h"
#include <mpi.h>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

// Start with an active box-only action. The strided atomic bias first becomes
// active on the following step, after the initial all-atom exchange.
bool delayedAtomicBias(MPI_Comm communicator, int rank, const char* asynchronous) {
  setenv("PLUMED_ASYNC_SHARE",asynchronous,1);
  PLMD::Plumed p;
  int precision=8, natoms=2, nlocal=1, index=rank;
  double timestep=0.001;
  p.cmd("setMPIComm",&communicator);
  p.cmd("setRealPrecision",&precision);
  p.cmd("setNatoms",&natoms);
  p.cmd("setTimestep",&timestep);
  p.cmd("setLogFile","/dev/null");
  p.cmd("init");
  p.cmd("readInputLine","v: VOLUME");
  p.cmd("readInputLine","PRINT ARG=v STRIDE=1 FILE=/dev/null");
  p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
  p.cmd("readInputLine","b: RESTRAINT ARG=d AT=0 KAPPA=1 STRIDE=2");
  p.cmd("setAtomsNlocal",&nlocal);
  p.cmd("setAtomsGatindex",&index);
  bool passed=true;
  for(int step=1; step<=4; ++step) {
    const double distance=1.0+step;
    std::array<double,3> positions{{1.0+rank*distance,2.0,3.0}};
    std::array<double,3> forces{{0.0,0.0,0.0}};
    std::array<double,9> box{{20,0,0,0,20,0,0,0,20}}, virial{};
    double mass=1.0, charge=0.0, bias=0.0;
    p.cmd("setStep",&step);
    p.cmd("setPositions",positions.data());
    p.cmd("setForces",forces.data());
    p.cmd("setMasses",&mass);
    p.cmd("setCharges",&charge);
    p.cmd("setBox",box.data());
    p.cmd("setVirial",virial.data());
    p.cmd("calc");
    p.cmd("getBias",&bias);
    const double expectedBias=step%2==0 ? 0.5*distance*distance : 0.0;
    const double expectedForce=step%2==0 ? (rank==0 ? 2.0 : -2.0)*distance : 0.0;
    const bool stepPassed=std::isfinite(bias) && std::abs(bias-expectedBias)<1e-12 &&
                          std::isfinite(forces[0]) && std::abs(forces[0]-expectedForce)<1e-12 &&
                          std::abs(forces[1])<1e-12 && std::abs(forces[2])<1e-12;
    if(!stepPassed) {
      std::cerr << "async=" << asynchronous << " rank=" << rank << " step=" << step
                << " bias=" << bias << " expected_bias=" << expectedBias
                << " force=" << forces[0] << " expected_force=" << expectedForce << '\n';
    }
    passed=passed && stepPassed;
  }
  return passed;
}

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  int rank=0,size=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(size!=2) {
    MPI_Finalize();
    return 2;
  }
  std::ofstream output;
  if(rank==0) {
    output.open("output");
  }
  int allPassed=1;
  for(const char* asynchronous : {
        "no","yes"
      }) {
    bool passed=false;
    try {
      passed=delayedAtomicBias(MPI_COMM_WORLD,rank,asynchronous);
    } catch(const std::exception& error) {
      std::cerr << error.what() << '\n';
    }
    int local=passed ? 1 : 0, global=0;
    MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    if(rank==0)
      output << (global ? "PASS" : "FAIL")
             << " firststep-inactive-domain async=" << asynchronous << '\n';
    allPassed=allPassed && global;
  }
  MPI_Finalize();
  return allPassed ? 0 : 1;
}
