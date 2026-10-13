#include "plumed/wrapper/Plumed.h"

#include <mpi.h>

#include <array>
#include <fstream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  bool localPassed=false;
  {
    try {
      PLMD::Plumed plumed;
      MPI_Comm communicator=MPI_COMM_WORLD;
      int precision=8;
      int natoms=32;
      double timestep=0.001;
      plumed.cmd("setMPIComm",&communicator);
      plumed.cmd("setRealPrecision",&precision);
      plumed.cmd("setNatoms",&natoms);
      plumed.cmd("setTimestep",&timestep);
      plumed.cmd("init");
      plumed.cmd("readInputLine",
                 "v: VORONOI_COORDINATION CENTERS=1-16 ASSIGNED=17-32 KAPPA=5 "
                 "REFERENCE=1 POWER=2 NOPBC");
      plumed.cmd("readInputLine","PRINT ARG=v FILE=/dev/null");

      std::vector<double> positions(32*3,0.0);
      for(unsigned i=0; i<16; ++i) {
        positions[3*i]=static_cast<double>(i);
        positions[3*(16+i)]=100.0+static_cast<double>(i);
      }
      positions[3*31]=15.0;
      std::vector<double> forces(positions.size(),0.0);
      std::vector<double> masses(32,1.0);
      std::array<double,9> box{{300,0,0,0,300,0,0,0,300}};
      std::array<double,9> virial{{0,0,0,0,0,0,0,0,0}};
      int step=11;
      plumed.cmd("setStep",&step);
      plumed.cmd("setPositions",positions.data());
      plumed.cmd("setMasses",masses.data());
      plumed.cmd("setBox",box.data());
      plumed.cmd("setForces",forces.data());
      plumed.cmd("setVirial",virial.data());
      plumed.cmd("calc");
    } catch(const PLMD::Plumed::Exception& error) {
      const std::string message=error.what();
      localPassed=
        message.find("category=distance-zero-or-near-zero")!=std::string::npos &&
        message.find("step=11")!=std::string::npos &&
        message.find("pair_index=255(zero-based)")!=std::string::npos &&
        message.find("center_atom=16 assigned_atom=32")!=std::string::npos &&
        message.find("rank=3")!=std::string::npos;
    }
  }
  int localValue=localPassed ? 1 : 0;
  int globalValue=0;
  MPI_Allreduce(&localValue,&globalValue,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if(rank==0) {
    std::ofstream output("output");
    output << (globalValue ? "PASS" : "FAIL")
           << " mpi-omp-deterministic-first-pair-no-deadlock\n";
  }
  MPI_Finalize();
  return globalValue ? 0 : 1;
}
