#include "plumed/wrapper/Plumed.h"
#include "mpi-test.h"
#if defined(__PLUMED_MPIWRAP)
#include "plumed/tools/Communicator.h"
#endif
#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

namespace PLMD {
namespace regtest_mpi {
int run(MPI_Comm world,int argc,char** argv) {
  int rank=0,size=0;
  MPI_Comm_rank(world,&rank);
  MPI_Comm_size(world,&size);
  if(size!=4) {
    MPI_Abort(world,2);
  }
  const int replica=rank/2,local=rank%2;
  const int negative=argc>1 ? std::atoi(argv[1]) : 0;
  MPI_Comm intra,inter;
  MPI_Comm_split(world,replica,local,&intra);
  MPI_Comm_split(world,local,replica,&inter);
#if defined(__PLUMED_MPIWRAP)
  {
    bool blocked=false;
    Communicator policy,derived;
    policy.Set_comm(inter);
    policy.setCommunicationCheck([&]() {
      if(blocked) {
        throw std::runtime_error("derived communication blocked");
      }
    });
    policy.Split(0,replica,derived);
    blocked=true;
    bool caught=false;
    try {
      derived.Barrier();
    } catch(const std::runtime_error& e) {
      caught=std::string(e.what())=="derived communication blocked";
    }
    if(!caught || derived.Get_size()!=2) {
      throw std::runtime_error("split communicator lost its communication policy");
    }
  }
#endif
  for(const std::string action : {
        "maxent","walkers","ensemble"
      }) {
    if(negative && action!="ensemble") {
      continue;
    }
    if(!negative && action=="ensemble") {
      continue;
    }
    Plumed p;
    p.cmd("GREX setMPIIntracomm",&intra);
    if(local==0) {
      p.cmd("GREX setMPIIntercomm",&inter);
    }
    p.cmd("GREX init");
    p.cmd("setMPIComm",&intra);
    p.cmd("setNatoms",2);
    p.cmd("setTimestep",0.001);
    p.cmd("setKbT",2.5);
    p.cmd("init");
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    if(action=="maxent") {
      p.cmd("readInputLine","m: MAXENT ARG=d TYPE=EQUAL AT=0 KAPPA=1 TAU=10 PACE=1 TEMP=300 FILE=LAGMULT");
    } else if(action=="walkers") {
      p.cmd("readInputLine","m: METAD ARG=d SIGMA=1 HEIGHT=1 PACE=1 WALKERS_MPI FILE=HILLS");
    } else {
      p.cmd("readInputLine","ensemble: ENSEMBLE ARG=d");
      p.cmd("readInputLine","r: RESTRAINT ARG=ensemble.d AT=0 KAPPA=2");
    }
    std::array<double,6> x{},f{};
    std::array<double,2> mass{{1,1}};
    std::array<double,9> virial{},box{{10,0,0,0,10,0,0,0,10}};
    x[3]=1+replica;
    auto bind=[&](int step) {
      f.fill(0);
      virial.fill(0);
      p.cmd("setStep",step);
      p.cmd("setPositions",x.data());
      p.cmd("setMasses",mass.data());
      p.cmd("setBox",box.data());
      p.cmd("setForces",f.data());
      p.cmd("setVirial",virial.data());
    };
    for(int step=0; step<2; ++step) {
      bind(step);
      p.cmd("calc");
      p.cmd("prepareCalc");
      p.cmd("performCalcNoForces");
      double own=0;
      p.cmd("getBias",&own);
      if(negative) {
        // Only one replica enters auxiliary mode: it must abort before the
        // unmatched inter-replica operation or the subsequent local broadcast.
        if(replica==0) {
          p.cmd("setEvaluationMode",negative);
          bind(step);
          p.cmd("prepareCalc");
          p.cmd(negative==1 ? "performCalcNoUpdate" : "performCalcNoForces");
        }
        MPI_Barrier(world);
        return 9;
      }
      for(int mode : {
            1,2
          }) {
        for(int repeat=0; repeat<=replica; ++repeat) {
          p.cmd("setEvaluationMode",mode);
          bind(step);
          p.cmd("prepareCalc");
          p.cmd(mode==1 ? "performCalcNoUpdate" : "performCalcNoForces");
          double bias=0;
          p.cmd("getBias",&bias);
          if(!std::isfinite(bias) || std::abs(bias-own)>1e-10) {
            throw std::runtime_error("local auxiliary call changed the current model");
          }
          p.cmd("setEvaluationMode",0);
        }
      }
      if(action=="maxent") {
        p.cmd("GREX savePositions");
        p.cmd("GREX setPartner",replica^1);
        p.cmd("GREX prepare");
        p.cmd("GREX calculate");
        double delta=0;
        p.cmd("GREX getLocalDeltaBias",&delta);
        const double expected=own*(double(2-replica)/double(1+replica)-1.0);
        if(!std::isfinite(delta) || std::abs(delta-expected)>1e-10) {
          throw std::runtime_error("MAXENT GREX mismatch rank="+std::to_string(rank)+" own="+std::to_string(own)+" delta="+std::to_string(delta)+" expected="+std::to_string(expected));
        }
      }
    }
  }
  MPI_Comm_free(&intra);
  MPI_Comm_free(&inter);
  if(rank==0) {
    std::ofstream("output") << "MAXENT GREX and unequal local auxiliary counts with MAXENT/WALKERS_MPI: PASS\n";
  }
  return 0;
}
}
}
int main(int argc,char** argv) {
  return PLMD::mpi_test::run(argc,argv,4,PLMD::regtest_mpi::run);
}
