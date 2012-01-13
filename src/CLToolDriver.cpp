#include "CLTool.h"
#include "CLToolRegister.h"
#include "Tools.h"
#include "Plumed.h"
#include <cstdio>
#include <string>
#include <vector>

using namespace std;

namespace PLMD {

/**
Class defining the driver

The driver is a tool to use plumed to process
an existing trajectory.
*/
class CLToolDriver:
public CLTool
{
public:
  int main(int argc,char**argv,FILE*in,FILE*out,PlumedCommunicator& pc);
  string description()const{
    return "analyze trajectories with plumed";
  }
};


PLUMED_REGISTER_CLTOOL(CLToolDriver,"driver")

int CLToolDriver::main(int argc,char**argv,FILE*in,FILE*out,PlumedCommunicator& pc){

// to avoid warnings:
 (void) in;

 string plumedFile("plumed.dat");
 string dumpforces("");
 string trajectoryFile("");
 double timestep(1.0);

 for(int i=1;i<argc;i++){
   string arg(argv[i]);
   if(arg.length()<1) continue;
   if(arg[0]=='-'){
     if(arg.find("--plumed=")==0){
       arg.erase(0,arg.find("=")+1);
       plumedFile=arg;
     } else if(arg.find("--timestep=")==0){
       arg.erase(0,arg.find("=")+1);
       Tools::convert(arg,timestep);
     } else if(arg.find("--dumpforces=")==0){
       arg.erase(0,arg.find("=")+1);
       dumpforces=arg;
       plumed_merror("not yet");
     } else plumed_merror("driver: unknown options "+arg);
   } else {
     if(trajectoryFile.length()==0) trajectoryFile=arg;
     else plumed_merror("driver: maximum one file at a time");
   }
 }
 if(trajectoryFile.length()==0) plumed_merror("please specify a trajectory");

  Plumed p;
  int checknatoms=0;
  int step=0;
  

  FILE* fp=fopen(trajectoryFile.c_str(),"r");
  
  std::string line;
  while(Tools::getline(fp,line)){

    int natoms;
    bool ok;
    Tools::convert(line,natoms);
    if(checknatoms==0){
      checknatoms=natoms;
      if(PlumedCommunicator::initialized()) p.cmd("setMPIComm",&pc.Get_comm());
      p.cmd("setNatoms",&natoms);
      p.cmd("setMDEngine","driver");
      p.cmd("setTimestep",&timestep);
      p.cmd("setPlumedDat",plumedFile.c_str());
      p.cmd("setLog",out);
      p.cmd("init");
    }
    plumed_massert(checknatoms==natoms,"number of atom changed");

    std::vector<double> coordinates(3*natoms,0.0);
    std::vector<double> forces(3*natoms,0.0);
    std::vector<double> masses(natoms,1.0);
    std::vector<double> cell(9,0.0);
    std::vector<double> virial(9,0.0);

    ok=Tools::getline(fp,line);
    plumed_massert(ok,"premature end of file");

    std::vector<std::string> words;
    words=Tools::getWords(line);
    plumed_massert(words.size()==3,"needed box in second line");
    for(unsigned i=0;i<3;i++) Tools::convert(words[0],cell[4*i]);
    for(unsigned i=0;i<natoms;i++){
      ok=Tools::getline(fp,line);
      plumed_massert(ok,"premature end of file");
      char dummy[1000];
// this is much faster than getWord
      std::sscanf(line.c_str(),"%s %lf %lf %lf",dummy,&coordinates[3*i],&coordinates[3*i+1],&coordinates[3*i+2]);
    }

   p.cmd("setForces",&forces[0]);
   p.cmd("setPositions",&coordinates[0]);
   p.cmd("setMasses",&masses[0]);
   p.cmd("setBox",&cell[0]);
   p.cmd("setStep",&step);
   p.cmd("calc");

    step++;
  }

  fclose(fp);
  
  return 0;
}



}
