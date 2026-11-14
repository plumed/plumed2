#include "plumed/wrapper/Plumed.h"
#include "mpi-test.h"
#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace PLMD {
namespace regtest_mpi {
int run(MPI_Comm world, int, char**) {
  int rank,size;
  MPI_Comm_rank(world,&rank);
  MPI_Comm_size(world,&size);
  if(size!=2) {
    MPI_Abort(world,2);
  }
  auto evaluate=[&](bool energy, bool direct, bool lambda=false, bool local=false, bool offStride=false, int inputPolicy=0) {
    Plumed p;
    MPI_Comm intra=MPI_COMM_SELF;
    p.cmd("GREX setMPIIntracomm",&intra);
    p.cmd("GREX setMPIIntercomm",&world);
    p.cmd("GREX init");
    p.cmd("setMPIComm",&intra);
    int natoms=2;
    double dt=.001;
    p.cmd("setNatoms",&natoms);
    p.cmd("setTimestep",&dt);
    p.cmd("init");
    p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    if(local) {
      if(inputPolicy==1) {
        p.cmd("readInputLine","lambda: PUT UNIT=number SHAPE=0 MUTABLE PERIODIC=NO NOEXCHANGE");
        p.cmd("readInputLine","l: COMBINE ARG=lambda PERIODIC=NO");
      } else {
        p.cmd("readInputLine","l: EXTRACV NAME=lambda");
      }
      p.cmd("readInputLine","c: CUSTOM ARG=l,d FUNC=x*y PERIODIC=NO");
      p.cmd("readInputLine","r: RESTRAINT ARG=c AT=0 KAPPA=2");
    }
    if(energy) {
      p.cmd("readInputLine","e: ENERGY");
      p.cmd("readInputLine",offStride ? "b: RESTRAINT ARG=e AT=0 KAPPA=2 STRIDE=2" : "b: RESTRAINT ARG=e AT=0 KAPPA=2");
    } else {
      const std::string input="b: METAD ARG=d SIGMA=1 HEIGHT="+std::to_string(.2*(rank+1))+" PACE=1 FILE=HILLS";
      p.cmd("readInputLine",input.c_str());
    }
    std::array<double,6> x{},f{};
    std::array<double,9> box{},virial{};
    std::array<double,2> masses{{1,1}};
    x[3]=rank+1;
    double own=(rank+1)*x[3]*x[3];
    const double otherDistance=2-rank;
    double foreign=(rank+1)*otherDistance*otherDistance;
    double parameter=rank+1;
    double parameterForce=0;
    for(int step=0; step<2; step++) {
      f.fill(0);
      virial.fill(0);
      p.cmd("setStep",&step);
      if(local) {
        parameterForce=0;
        p.cmd(inputPolicy==0 ? "setValue lambda NOEXCHANGE" : "setExtraCV lambda",&parameter);
        p.cmd("setExtraCVForce lambda",&parameterForce);
      }
      p.cmd("setPositions",x.data());
      p.cmd("setMasses",masses.data());
      p.cmd("setBox",box.data());
      p.cmd("setForces",f.data());
      p.cmd("setVirial",virial.data());
      p.cmd("prepareCalc");
      if(energy) {
        p.cmd("setEnergy",&own);
      }
      p.cmd("performCalc");
    }
    double delta;
    if(direct) {
      const auto savedForces=f;
      const auto savedVirial=virial;
      int exchange=1, step=1;
      p.cmd("setStep",&step);
      p.cmd("setExchangeStep",&exchange);
      p.cmd("prepareCalc");
      p.cmd("performCalcNoForces");
      double localBias;
      p.cmd("getBias",&localBias);
      x[3]=otherDistance;
      p.cmd("setStep",&step);
      p.cmd("setPositions",x.data());
      if(energy) {
        p.cmd("setEnergy",&foreign);
      }
      p.cmd("prepareCalc");
      p.cmd("performCalcNoForces");
      double foreignBias;
      p.cmd("getBias",&foreignBias);
      exchange=0;
      p.cmd("setExchangeStep",&exchange);
      if(f!=savedForces || virial!=savedVirial) {
        throw std::runtime_error("Energy-only evaluation changed forces");
      }
      delta=foreignBias-localBias+(energy ? foreign-own : 0);
    } else {
      p.cmd("GREX savePositions");
      if(energy) {
        p.cmd("GREX cacheLocalUNow",&own);
        p.cmd("GREX cacheLocalUSwap",&foreign);
      }
      int partner=rank^1;
      p.cmd("GREX setPartner",&partner);
      p.cmd("GREX calculate");
      p.cmd("GREX getLocalDeltaBias",&delta);
      if(energy && lambda) {
        delta-=foreign-own;
      }
    }
    const double gaussian=(std::exp(-.5)-std::exp(-6.25))/(1-std::exp(-6.25));
    double expected=energy ? foreign*foreign-own*own+(lambda ? 0 : foreign-own) : .2*(rank+1)*(gaussian-1);
    if(local) {
      const double foreignParameter=inputPolicy==2 ? 2-rank : parameter;
      expected+=foreignParameter*foreignParameter*otherDistance*otherDistance-parameter*parameter*x[3]*x[3];
    }
    if(std::abs(delta-expected)>1e-10) {
      throw std::runtime_error("GREX energy mismatch");
    }
  };
  evaluate(true,false);
  evaluate(false,false);
  evaluate(true,true);
  evaluate(false,true);
  evaluate(true,false,true);
  evaluate(false,false,false,true);
  evaluate(false,false,false,true,false,1); // PUT metadata retains lambda locally.
  evaluate(false,false,false,true,false,2); // Default extra CVs follow the coordinates.
  evaluate(true,false,true,false,true);
  MPI_Barrier(world);
  if(rank==0) {
    std::ofstream("output")<<"ENERGY cross-Hamiltonian evaluation: pass\nPost-update bias consistency: pass\nDirect energy-only evaluation: pass\nLambda ENERGY without physical double counting: pass\nLocal lambda retained on foreign coordinates: pass\nOff-stride exchange bias evaluation: pass\n";
  }
  return 0;
}
}
}
int main(int argc,char**argv) {
  return PLMD::mpi_test::run(argc,argv,2,PLMD::regtest_mpi::run);
}
