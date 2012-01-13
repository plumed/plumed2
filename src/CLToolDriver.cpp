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
 double timestep(0.001);
 unsigned stride(1);
 bool printhelp=false;

// Start parsing options
  string prefix("");
  string a("");
  for(int i=1;i<argc;i++){
    a=prefix+argv[i];
    if(a.length()==0) continue;
    if(a=="-h" || a=="--help"){
      printhelp=true;
      break;
    }
    if(a.find("--plumed=")==0){
      a.erase(0,a.find("=")+1);
      plumedFile=a;
      prefix="";
    } else if(a=="--plumed"){
      prefix="--plumed=";
    } else if(a.find("--timestep=")==0){
      a.erase(0,a.find("=")+1);
      Tools::convert(a,timestep);
      prefix="";
    } else if(a=="--timestep"){
      prefix="--timestep=";
    } else if(a.find("--stride=")==0){
      a.erase(0,a.find("=")+1);
      Tools::convert(a,stride);
      prefix="";
    } else if(a=="--stride"){
      prefix="--stride=";
    } else if(a.find("--dumpforces=")==0){
      a.erase(0,a.find("=")+1);
      dumpforces=a;
      plumed_merror("not yet");
    } else if(a[0]=='-') {
      string msg="ERROR: Unknown option " +a;
      fprintf(stderr,"%s\n",msg.c_str());
      return 1;
    } else if(trajectoryFile.length()==0){
      trajectoryFile=a;
    } else {
      string msg="ERROR: maximum one file at a time";
      fprintf(stderr,"%s\n",msg.c_str());
      return 1;
    }
  }

  if(printhelp){
    fprintf(out,"%s",
 "Usage: driver [options] trajectory.xyz\n"
 "Options:\n"
 "  [--help|-h]     : prints this help\n"
 "  [--plumed FILE] : plumed script file (default: plumed.dat)\n"
 "  [--timestep TS] : timestep (default: 0.001) in picoseconds\n"
 "  [--stride ST]   : stride between frames (default: 1)\n"
);
    return 0;
  }

  if(trajectoryFile.length()==0){
    string msg="ERROR: please specify a trajectory";
    fprintf(stderr,"%s\n",msg.c_str());
    return 1;
  }

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

    step+=stride;
  }

  fclose(fp);
  
  return 0;
}



}
