#include <fstream>

#include "plumed/wrapper/Plumed.h"
#include "plumed/tools/Communicator.h"

constexpr unsigned nat=10;
//it is a struct because we don't need getter/setters (for now)
struct plumedThrowChecker{
  unsigned natoms=nat;
  std::vector<double> positions=[](unsigned natoms){
    std::vector<double> toret(natoms*3,0.0);
    for(unsigned i=0; i<3*natoms; i++) toret[i]=i;
    return toret;
  }(nat);
  std::vector<double> masses{nat,1.0};
  std::vector<double> forces{3*nat,0.0};
  std::vector<double> box{9,0.0};
  std::vector<double> virial{9,0.0};
  ///This initializes a new plumed to test each throw in the cleanliest way possible
  int checkThrow (std::ostream& out,std::string name, std::string cmd, std::string expectedMessage){
    //GIVEN an initialized Plumed interface
    PLMD::Plumed plumed;
    
    plumed.cmd("setNatoms",&natoms);
    int step=0;
    plumed.cmd("setLogFile","test.log");
    plumed.cmd("init");
    plumed.cmd("setStep",&step);
    plumed.cmd("setPositions",positions.data());
    plumed.cmd("setBox",box.data());
    plumed.cmd("setForces",forces.data());
    plumed.cmd("setVirial",virial.data());
    plumed.cmd("setMasses",masses.data());

    plumed.cmd("readInputLine","d: DISTANCE ATOMS=1,2");
    plumed.cmd("readInputLine","d1: DISTANCE ATOMS={1 2}");
    ///TODO: expand with a "readInputLines" to give the possibility to test in more situations
    try {
      //WHEN the user ask for the given input
      plumed.cmd("readInputLine",cmd.c_str());
      //THEN plumed should gracefully exit with a clear error message
    } catch(PLMD::Plumed::ExceptionError &e) { //correct throw, we are happy
      std::string exceptionText{e.what()};
      out << name << " : ";
      if (exceptionText.find(expectedMessage) != std::string::npos) {
        out << "Exception thrown, correct message\n";
        return 0;
      }

      out << "Exception thrown, wrong message: "
                << e.what() ;
      out << "\tExpected message should contain: \""
                << expectedMessage << "\"\n";
      return 1;
    }
    out  << "Exception not thrown\n";
    return 1;
  }
};

int main(int, char**) {
  std::ofstream out("output");
  plumedThrowChecker ptc;
  //When the user aks for a WALKERS_MPI in the METAD action,
  //if MPI is not installed then the user must be informed
  //if MPI is installed then the comunications must be already set up
  //WHEN PLUMED is not compiled with MPI or the MPI routines are not initialized
  std::string expectedMessage="WALKERS_MPI flag requires MPI compilation";
  if (PLMD::Communicator::plumedHasMPI()) {
    expectedMessage="WALKERS_MPI needs the communicator correctly initialized";
  }
  ptc.checkThrow(out,"METAD WALKER_MPI",
    "METAD ARG=d,d1 SIGMA=0.1,0.2 HEIGHT=0.1 PACE=2 RESTART=YES WALKERS_MPI",
    expectedMessage);

  //add other throw messages checks
}
