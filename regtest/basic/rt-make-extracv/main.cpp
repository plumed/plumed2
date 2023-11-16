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

  plumed->cmd("setNatoms",natoms);
  plumed->cmd("setLogFile","test.log");
  plumed->cmd("init");
  plumed->cmd("readInputLine","d: DISTANCE ATOMS=1,2");
  plumed->cmd("readInputLine","e: EXTRACV NAME=extra");
  plumed->cmd("readInputLine","e2: EXTRACV NAME=extra2");
  plumed->cmd("readInputLine","e3: EXTRACV NAME=extra3");
  plumed->cmd("readInputLine","PRINT ARG=e FILE=COLVARX");
  plumed->cmd("readInputLine","RESTRAINT ARG=e AT=0 KAPPA=1");
  plumed->cmd("readInputLine","RESTRAINT ARG=d AT=0 KAPPA=1");
  plumed->cmd("readInputLine","RESTRAINT ARG=e2 AT=0 KAPPA=-1");
  plumed->cmd("readInputLine","PRINT ARG=e3 FILE=extra3 STRIDE=3");

  std::ofstream ofs("output");

  for(int step=0;step<10;step++){
    double extracv=step;
    double extracvf=0.0;
    double extracv2=step;
    double extracvf2=-step;
    double extracv3=2*step;
    double extracvf3=0.0;
    plumed->cmd("setStep",step);
    plumed->cmd("setPositions",&positions[0],3*natoms);
    plumed->cmd("setBox",&box[0],9);
    plumed->cmd("setForces",&forces[0],3*natoms);
    plumed->cmd("setVirial",&virial[0],9);
    plumed->cmd("setExtraCV extra",&extracv,1);
    plumed->cmd("setExtraCVForce extra",&extracvf,1);
    plumed->cmd("setExtraCV extra2",&extracv2,1);
    plumed->cmd("setExtraCVForce extra2",&extracvf2,1);
    plumed->cmd("setExtraCV extra3",&extracv3,1);
    plumed->cmd("setExtraCVForce extra3",&extracvf3,1);
    plumed->cmd("setMasses",&masses[0],natoms);
// first compute using modified positions:
    positions[0]=0.5;
    extracv2=100;
    for(auto & f:forces) f=0.0;
    extracvf=0.0;
    extracvf2=0.0;
    plumed->cmd("prepareCalc");
    int isExtraCV3Needed=0;
    plumed->cmd("isExtraCVNeeded extra3",&isExtraCV3Needed);
    ofs<<"extracv3_needed: "<<isExtraCV3Needed<<"\n";
    if(!isExtraCV3Needed) extracv3=0.0;
    plumed->cmd("performCalcNoUpdate");
    double bias;
    plumed->cmd("getBias",&bias);
    ofs<<"bias_pre: "<<bias<<"\n";
    ofs<<"extracvf_pre: "<<extracvf<<"\n";
    ofs<<"extracvf2_pre: "<<extracvf2<<"\n";
    ofs<<"f_pre:";
    for(auto & f:forces) ofs<<" "<<f;
    ofs<<"\n";
// thent compute using another modified positions, without updating forces
    positions[0]=0.5;
    extracv2=100;
    for(auto & f:forces) f=0.0;
    extracvf=0.0;
    extracvf2=0.0;
    plumed->cmd("prepareCalc");
    plumed->cmd("performCalcNoForces");
    plumed->cmd("getBias",&bias);
    ofs<<"bias_pre: "<<bias<<"\n";
    ofs<<"extracvf_pre: "<<extracvf<<"\n";
    ofs<<"extracvf2_pre: "<<extracvf2<<"\n";
    ofs<<"f_pre:";
    for(auto & f:forces) ofs<<" "<<f;
    ofs<<"\n";
// then compute using regular positions:
    positions[0]=0;
    extracv2=0;
    for(auto & f:forces) f=0.0;
    extracvf=0.0;
    extracvf2=0.0;
    plumed->cmd("prepareCalc");
    plumed->cmd("performCalcNoUpdate");
    plumed->cmd("getBias",&bias,1);
    ofs<<"bias: "<<bias<<"\n";
    ofs<<"extracvf: "<<extracvf<<"\n";
    ofs<<"extracvf2: "<<extracvf2<<"\n";
    ofs<<"f:";
    for(auto & f:forces) ofs<<" "<<f;
    ofs<<"\n";
    plumed->cmd("update");
  }

  delete plumed;
  return 0;
}
