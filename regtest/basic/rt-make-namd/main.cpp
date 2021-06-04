#include <vector>
#include <fstream>
#include "plumed/wrapper/Plumed.h"
#include <iostream>

using namespace PLMD;

int main(){
  Plumed* plumed=new Plumed;

  int natoms=100;

  std::vector<double> all_positions(3*natoms,0.0);
  for(unsigned i=0;i<natoms;i++) all_positions[i]=i;
  std::vector<double> all_masses(natoms,1.0);
  std::vector<double> all_charges(natoms,1.0);
  std::vector<double> all_forces(3*natoms,0.0);
  std::vector<double> box(9,0.0);
  std::vector<double> virial(9,0.0);

  plumed->cmd("setNatoms",&natoms);
  plumed->cmd("setLogFile","test.log");
  double dt=0.001;
  plumed->cmd("setTimestep",&dt);
  plumed->cmd("setNoVirial");

  std::string file="plumed.dat";
  plumed->cmd("setPlumedDat",file.c_str());

  plumed->cmd("init");
  std::ofstream ofs("output");

  std::vector<int> index;
  std::vector<double> masses;
  std::vector<double> forces;
  std::vector<double> positions;
  std::vector<double> charges;

  bool first=true;

  for(int step=0;step<10;step++){
  for(unsigned i=0;i<3*natoms;i++) all_positions[i]=i+step;
    const int* p=nullptr;
    int n=0;
    plumed->cmd("setStep",&step);
    plumed->cmd("prepareDependencies");
    plumed->cmd("createFullList",&n);
//std::cerr<<"n= "<<n<<"\n";
    plumed->cmd("getFullList",&p);
    bool redo=(index.size()!=n);
    if(first) redo=true;
    first=false;
    if(!redo) for(int i=0;i<n;i++) if(index[i]!=p[i]) { redo=true; break;};
    if(redo){
      index.resize(n);
      masses.resize(n);
      for(int i=0;i<n;i++){
        masses[i]=all_masses[p[i]];  
        index[i]=p[i];
      };
      positions.resize(3*n);
      forces.resize(3*n);
      charges.resize(n);
      
      plumed->cmd("setAtomsNlocal",&n);
      plumed->cmd("setAtomsGatindex",(index.empty()?nullptr:&index[0]));
// for(unsigned i=0;i<n;i++) std::cerr<< "I "<<index[i] <<"\n";
    }
    plumed->cmd("clearFullList");

    for(int i=0;i<index.size();i++){
      positions[3*i+0]=all_positions[3*index[i]+0];
      positions[3*i+1]=all_positions[3*index[i]+1];
      positions[3*i+2]=all_positions[3*index[i]+2];
      masses[i]=all_masses[index[i]];
      charges[i]=all_charges[index[i]];
    };

    plumed->cmd("setBox",&box[0]);

    for(int i=0;i<forces.size();i++) forces[i]=0.0;

// std::cerr<<"mass "<<&masses[0]<<" "<<(&masses[0]==NULL)<<"\n";
    plumed->cmd("setMasses",(masses.empty()?nullptr:&masses[0]));
    plumed->cmd("setCharges",(charges.empty()?nullptr:&charges[0]));
    plumed->cmd("setPositions",(positions.empty()?nullptr:&positions[0]));
    plumed->cmd("setForces",(forces.empty()?nullptr:&forces[0]));

    plumed->cmd("shareData");
    plumed->cmd("performCalc");

    std::vector<double> all_forces(3*natoms,0.0);

    for(unsigned i=0;i<all_forces.size();i++) all_forces[i]=0.0;
    for(int i=0;i<index.size();i++){
      all_forces[3*index[i]+0]=forces[3*i+0];
      all_forces[3*index[i]+1]=forces[3*i+1];
      all_forces[3*index[i]+2]=forces[3*i+2];
    }
  }

  delete plumed;
  return 0;
}
