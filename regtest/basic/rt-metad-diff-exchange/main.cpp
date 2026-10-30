#include "plumed/wrapper/Plumed.h"
#include <mpi.h>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

// distance followed by configuration c at step t
static double distance(int c, int t) {
  return 1.0+2.0*c+0.05*std::sin(0.3*t+c)+0.03*std::sin(1.7*t);
}

static std::string slurp(const std::string& name) {
  std::ifstream f(name);
  std::stringstream s;
  s<<f.rdbuf();
  return s.str();
}

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(size!=4) {
    MPI_Abort(MPI_COMM_WORLD,2);
  }
  MPI_Comm intra=MPI_COMM_SELF;
  MPI_Comm inter=MPI_COMM_WORLD;
  int partner=rank^1;
  {
    PLMD::Plumed p;
    const int natoms=2;
    const double timestep=0.002;
    p.cmd("GREX setMPIIntracomm",&intra);
    p.cmd("GREX setMPIIntercomm",&inter);
    p.cmd("GREX init");
    p.cmd("setMPIComm",&intra);
    p.cmd("setNatoms",&natoms);
    p.cmd("setTimestep",&timestep);
    const std::string logfile="test."+std::to_string(rank)+".log";
    p.cmd("setLogFile",logfile.c_str());
    p.cmd("init");
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    p.cmd("readInputLine","m: METAD ARG=d ADAPTIVE=DIFF SIGMA=20 HEIGHT=0.0 PACE=10 FILE=HILLS FMT=%10.5f");
    std::array<double,6> x{}, f{};
    std::array<double,9> box{}, virial{};
    std::array<double,2> masses{{1,1}};
    int conf=rank%2; // configuration followed by this replica
    for(int step=0; step<300; step++) {
      const bool exchange=(step==100) || (step==200 && rank<2);
      x[3]=distance(conf,step);
      f.fill(0);
      virial.fill(0);
      p.cmd("setStep",&step);
      p.cmd("setPositions",x.data());
      p.cmd("setMasses",masses.data());
      p.cmd("setBox",box.data());
      p.cmd("setForces",f.data());
      p.cmd("setVirial",virial.data());
      p.cmd("prepareCalc");
      if(exchange) {
        p.cmd("GREX savePositions");
      }
      p.cmd("performCalc");
      if(exchange) {
        p.cmd("GREX setPartner",&partner);
        p.cmd("GREX calculate");
        // the exchange at step 100 is accepted, the one at step 200 is rejected
        if(step==100) {
          conf=1-conf;
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0) {
    std::ofstream out("output");
    out<<"rejected exchange has no effect on replica 0: "<<(slurp("HILLS.0")==slurp("HILLS.2") ? "yes" : "no")<<"\n";
    out<<"rejected exchange has no effect on replica 1: "<<(slurp("HILLS.1")==slurp("HILLS.3") ? "yes" : "no")<<"\n";
  }
  MPI_Finalize();
  return 0;
}
