#include "plumed/wrapper/Plumed.h"
#include "mpi-test.h"
#include <cmath>
#include <fstream>
#include <string>

namespace PLMD {
namespace regtest_mpi {
int run(MPI_Comm world, int, char**) {
  int rank=0,size=0;
  MPI_Comm_rank(world,&rank);
  MPI_Comm_size(world,&size);
  bool passed=true;
  for(int beads : {1,2}) {
    if(size%beads) { return 2; }
    const int bead=rank/(size/beads), local=rank%(size/beads);
    MPI_Comm intra,inter;
    MPI_Comm_split(world,bead,local,&intra);
    MPI_Comm_split(world,local,bead,&inter);
    for(bool weighted : {false,true}) {
      for(bool central : {false,true}) {
        for(bool powered : {false,true}) {
          Plumed p;
          p.cmd("GREX setMPIIntracomm",&intra);
          if(local==0) { p.cmd("GREX setMPIIntercomm",&inter); }
          p.cmd("GREX init");
          p.cmd("setMPIComm",&intra);
          const int natoms=0;
          const double dt=1, kbt=1;
          p.cmd("setNatoms",&natoms);
          p.cmd("setTimestep",&dt);
          p.cmd("setKbT",&kbt);
          p.cmd("init");
          const std::string c="c: CONSTANT VALUE="+std::to_string(1+2*bead);
          p.cmd("readInputLine",c.c_str());
          const std::string e="energy: CONSTANT VALUE="+std::to_string(bead*std::log(3.0));
          p.cmd("readInputLine",e.c_str());
          std::string input="a: ENSEMBLE ARG=c";
          if(weighted) { input+=",energy REWEIGHT"; }
          input+=" MOMENT=2";
          if(central) { input+=" CENTRAL"; }
          if(powered) { input+=" POWER=2"; }
          p.cmd("readInputLine",input.c_str());
          p.cmd("readInputLine","b: BIASVALUE ARG=a.c,a.c_m");
          const int step=0;
          p.cmd("setStep",&step);
          p.cmd("calc");
          double bias=0;
          p.cmd("getBias",&bias);
          const double w=beads==1 ? 0 : (weighted ? 0.75 : 0.5);
          double mean=1+2*w, moment=central ? 4*w*(1-w) : 1+8*w;
          if(powered) { mean*=mean; moment*=moment; }
          passed=passed && std::abs(bias-mean-moment)<1e-5;
        }
      }
    }
    MPI_Comm_free(&intra);
    MPI_Comm_free(&inter);
  }
  int ok=passed,all=0;
  MPI_Allreduce(&ok,&all,1,MPI_INT,MPI_MIN,world);
  if(rank==0) { std::ofstream("output") << (all ? "PASS" : "FAIL") << "\n"; }
  return all ? 0 : 1;
}
}
}
int main(int argc, char** argv) {
  return PLMD::mpi_test::run(argc,argv,4,PLMD::regtest_mpi::run);
}
