#include "plumed/wrapper/Plumed.h"
#include "plumed/function/Function.h"
#include "plumed/core/ActionRegister.h"
#include "plumed/tools/Communicator.h"
#include <mpi.h>
#include <limits>
#include <array>
#include <cstdlib>
#include <fstream>
#include <string>

// The same action graph is used on all ranks. Only the synthetic derivative
// is rank-local, as can happen during distributed force calculations.
namespace PLMD {
class RankLocalTraceTest : public function::Function {
  bool inject;
public:
  static void registerKeywords(Keywords& keys) {
    function::Function::registerKeywords(keys);
    keys.addFlag("INJECT",false,"inject a non-finite derivative on rank three");
  }
  explicit RankLocalTraceTest(const ActionOptions& options)
    : Action(options), function::Function(options), inject(false) {
    parseFlag("INJECT",inject);
    checkRead();
    addValueWithDerivatives();
    setNotPeriodic();
  }
  void calculate() override {
    setValue(getArgument(0));
    setDerivative(0,inject && comm.Get_rank()==3 ?
                  std::numeric_limits<double>::infinity() : 1.0);
  }
};
PLUMED_REGISTER_ACTION(RankLocalTraceTest,"RANK_LOCAL_TRACE_TEST")
}

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  setenv("PLUMED_NONFINITE_ACTION_TRACE","c",1);
  std::ofstream output;
  if(rank==0) {
    output.open("output");
  }
  bool passed=true;
  for(int singular=1; singular>=0; --singular) {
    std::string message;
    try {
      PLMD::Plumed p;
      MPI_Comm communicator=MPI_COMM_WORLD;
      int natoms=2;
      double timestep=0.001;
      p.cmd("setMPIComm",&communicator);
      p.cmd("setNatoms",&natoms);
      p.cmd("setTimestep",&timestep);
      p.cmd("init");
      p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
      const std::string input=singular ? "c: RANK_LOCAL_TRACE_TEST ARG=d INJECT" :
                              "c: RANK_LOCAL_TRACE_TEST ARG=d";
      p.cmd("readInputLine",input.c_str());
      p.cmd("readInputLine","r: RESTRAINT ARG=c AT=1 KAPPA=1");
      std::array<double,6> positions{{0,0,0,1,0,0}},forces{};
      std::array<double,2> masses{{1,1}};
      std::array<double,9> box{{10,0,0,0,10,0,0,0,10}},virial{};
      int step=17;
      p.cmd("setStep",&step);
      p.cmd("setPositions",positions.data());
      p.cmd("setMasses",masses.data());
      p.cmd("setBox",box.data());
      p.cmd("setForces",forces.data());
      p.cmd("setVirial",virial.data());
      p.cmd("calc");
    } catch(const PLMD::Plumed::Exception& error) {
      message=error.what();
    }
    int success=singular ? (message.find("rank=3")!=std::string::npos &&
                            message.find("field=derivative")!=std::string::npos)
                : message.empty();
    MPI_Allreduce(MPI_IN_PLACE,&success,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    passed=passed && success;
    if(rank==0)
      output << (success ? "PASS " : "FAIL ")
             << (singular ? "rank-local-error-propagation" : "finite-rank-local-values") << "\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
