#include "plumed/wrapper/Plumed.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

// Two atoms may belong to either rank or be shared between the ranks.
// Keep dummy MD storage on the empty rank; PLUMED must handle empty messages.
int run(MPI_Comm communicator, int rank, const char* asynchronous, int owner,
        std::ofstream& output) {
  if(setenv("PLUMED_ASYNC_SHARE",asynchronous,1)!=0) {
    MPI_Abort(communicator,2);
  }
  PLMD::Plumed p;
  int precision=sizeof(double), natoms=2;
  const int nlocal=owner==2 ? 1 : (owner==rank ? 2 : 0);
  int indices[2]= {owner==2 ? rank : 0,1};
  double timestep=0.001;
  p.cmd("setMPIComm",&communicator);
  p.cmd("setRealPrecision",&precision);
  p.cmd("setNatoms",&natoms);
  p.cmd("setTimestep",&timestep);
  p.cmd("setLogFile","/dev/null");
  p.cmd("init");
  p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
  p.cmd("readInputLine","b: RESTRAINT ARG=d AT=1 KAPPA=1");
  p.cmd("setAtomsNlocal",&nlocal);
  p.cmd("setAtomsGatindex",indices);
  int failures=0;
  for(int step=0; step<2; ++step) {
    double positions[2][3]= {{2.0*indices[0],0,0},{2,0,0}};
    double forces[2][3]= {};
    double masses[2]= {1,1}, charges[2]= {0,0};
    double box[3][3]= {{10,0,0},{0,10,0},{0,0,10}}, virial[3][3]= {};
    double bias=0;
    p.cmd("setStep",&step);
    p.cmd("setPositions",&positions[0][0], {2,3});
    p.cmd("setForces",&forces[0][0], {2,3});
    p.cmd("setMasses",masses);
    p.cmd("setCharges",charges);
    p.cmd("setBox",&box[0][0], {3,3});
    p.cmd("setVirial",&virial[0][0], {3,3});
    p.cmd("calc");
    p.cmd("getBias",&bias);
    bool passed=std::isfinite(bias) && std::abs(bias-0.5)<1e-12;
    double globalForces[2][3]= {};
    for(int atom=0; atom<nlocal; ++atom) {
      for(int axis=0; axis<3; ++axis) {
        globalForces[indices[atom]][axis]=forces[atom][axis];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,globalForces,6,MPI_DOUBLE,MPI_SUM,communicator);
    for(int atom=0; atom<2; ++atom) {
      for(int axis=0; axis<3; ++axis) {
        const double expected=axis==0 ? (atom==0 ? 1.0 : -1.0) : 0.0;
        passed=passed && std::isfinite(globalForces[atom][axis]) &&
               std::abs(globalForces[atom][axis]-expected)<1e-12;
      }
    }
    int local=passed ? 1 : 0, global=0;
    MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,communicator);
    if(rank==0) {
      output << (global ? "PASS" : "FAIL") << " empty-rank async=" << asynchronous
             << " owner=" << owner << " step=" << step << '\n';
      output.flush();
    }
    if(!global) {
      ++failures;
      std::cerr << "rank=" << rank << " async=" << asynchronous << " owner=" << owner
                << " step=" << step << " bias=" << bias
                << " forces=" << globalForces[0][0] << ',' << globalForces[1][0] << '\n';
    }
  }
  return failures;
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
    if(!output) {
      MPI_Abort(MPI_COMM_WORLD,2);
    }
  }
  int failures=0;
  try {
    for(const char* asynchronous : {
          "no","yes"
        }) {
      for(int owner=0; owner<3; ++owner) {
        failures+=run(MPI_COMM_WORLD,rank,asynchronous,owner,output);
      }
    }
  } catch(const std::exception& error) {
    std::cerr << "rank=" << rank << " " << error.what() << '\n';
    MPI_Abort(MPI_COMM_WORLD,3);
  }
  MPI_Finalize();
  return failures ? 1 : 0;
}
