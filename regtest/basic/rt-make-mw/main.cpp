#include "plumed/wrapper/Plumed.h"
#include <vector>
#include <sstream>

using namespace PLMD;

void go(Plumed p,int natoms,unsigned iw,unsigned is){
  std::vector<double> positions(3*natoms,0.0);
  for(unsigned i=0;i<natoms;i++) positions[i]=i+iw+is;
  std::vector<double> masses(natoms,1.0);
  std::vector<double> forces(3*natoms,0.0);
  std::vector<double> box(9,0.0);
  std::vector<double> virial(9,0.0);

  p.cmd("setStep",&is);
  p.cmd("setPositions",&positions[0]);
  p.cmd("setBox",&box[0]);
  p.cmd("setForces",&forces[0]);
  p.cmd("setVirial",&virial[0]);
  p.cmd("setMasses",&masses[0]);
  p.cmd("calc");
}

int main(){
  std::vector<Plumed> p;

  unsigned nwalkers=3;
  unsigned nsteps=10;

  p.resize(nwalkers);

  int natoms=10;

  for(unsigned iw=0;iw<nwalkers;iw++){

    p[iw].cmd("setNatoms",&natoms);

    std::ostringstream iwss;
    iwss<<iw;
    std::string file;
    file="test." + iwss.str() + ".log";
    p[iw].cmd("setLogFile",file.c_str());
    file="plumed." + iwss.str() + ".dat";
    p[iw].cmd("setPlumedDat",file.c_str());
    p[iw].cmd("init");
  }

// half steps for each walker
  for(unsigned iw=0;iw<nwalkers;iw++) for(unsigned is=0;is<nsteps/2;is++) go(p[iw],natoms,iw,is);

// other half steps for each walker
  for(unsigned iw=0;iw<nwalkers;iw++) for(unsigned is=nsteps/2;is<nsteps;is++) go(p[iw],natoms,iw,is);

  return 0;
}
