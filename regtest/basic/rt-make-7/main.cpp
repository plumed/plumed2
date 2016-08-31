#include "plumed/wrapper/Plumed.h"
#include <vector>
#include <fstream>

using namespace PLMD;

int main(){
  Plumed* plumed=new Plumed;

  int natoms=10;

  std::vector<double> positions(3*natoms,0.0);
  for(unsigned i=0;i<natoms;i++) positions[i]=i;
  std::vector<double> masses(natoms,1.0);
  std::vector<double> forces(3*natoms,0.0);
  std::vector<double> box(9,0.0);
  std::vector<double> virial(9,0.0);

  plumed->cmd("setNatoms",&natoms);
  plumed->cmd("setLogFile","test.log");
  plumed->cmd("init");
  plumed->cmd("readInputLine","d: DISTANCE ATOMS=1,2");
  plumed->cmd("readInputLine","d1: DISTANCE ATOMS={1 2}"); // check if braces are parsed correctly
  plumed->cmd("readInputLine","PRINT ARG=d,d1 FILE=COLVAR");
  plumed->cmd("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1");
  plumed->cmd("readInputLine","METAD ARG=d PACE=1 SIGMA=1 HEIGHT=0 FILE=H1");
  plumed->cmd("readInputLine","METAD ARG=d PACE=2 SIGMA=1 HEIGHT=0 FILE=H2");

  std::ofstream ofs("output");

  for(int step=0;step<10;step++){
    plumed->cmd("setStep",&step);
    plumed->cmd("setPositions",&positions[0]);
    plumed->cmd("setBox",&box[0]);
    plumed->cmd("setForces",&forces[0]);
    plumed->cmd("setVirial",&virial[0]);
    plumed->cmd("setMasses",&masses[0]);
// first compute using modified positions:
    positions[0]=0.5;
    plumed->cmd("prepareCalc");
    plumed->cmd("performCalcNoUpdate");
    positions[0]=0;
    double bias;
    plumed->cmd("getBias",&bias);
    ofs<<bias<<"\n";
// first compute using regular positions:
    plumed->cmd("prepareCalc");
    plumed->cmd("performCalcNoUpdate");
    plumed->cmd("getBias",&bias);
    ofs<<bias<<"\n";
// hills should only be added at regular positions:
    plumed->cmd("update");
  }

  delete plumed;
  return 0;
}
