#include "CLTool.h"
#include "CLToolRegister.h"
#include "Tools.h"
#include "Plumed.h"
#include "PlumedCommunicator.h"
#include <cstdio>
#include <string>
#include <vector>

using namespace std;

namespace PLMD {

//+PLUMEDOC TOOLS driver
/**
driver is a tool that allows one to to use plumed to post-process
an existing trajectory.
*/
//+ENDPLUMEDOC

template<typename real>
class CLToolDriver:
public CLTool
{
public:
  int main(int argc,char**argv,FILE*in,FILE*out,PlumedCommunicator& pc);
  string description()const;
};

template<typename real>
string CLToolDriver<real>::description()const{ return "analyze trajectories with plumed"; }

template<>
string CLToolDriver<float>::description()const{ return "analyze trajectories with plumed (single precision version)"; }


template<typename real>
int CLToolDriver<real>::main(int argc,char**argv,FILE*in,FILE*out,PlumedCommunicator& pc){

// to avoid warnings:
 (void) in;

 string plumedFile("plumed.dat");
 string dumpforces("");
 string dumpforcesFmt("%f");
 string trajectoryFile("");
 real timestep(real(0.001));
 unsigned stride(1);
 bool printhelp=false;
 bool printhelpdebug=false;

// Start parsing options
  string prefix("");
  string a("");
  for(int i=1;i<argc;i++){
    a=prefix+argv[i];
    if(a.length()==0) continue;
    if(a=="-h" || a=="--help"){
      printhelp=true;
      break;
    } else if(a=="--help-debug"){
      printhelp=true;
      printhelpdebug=true;
      break;
    } else if(a=="--debug-float"){
      if(sizeof(real)!=sizeof(float)){
        CLTool* cl=new CLToolDriver<float>;
        int ret=cl->main(argc,argv,in,out,pc);
        delete cl;
        return ret;
      }
    } else if(a.find("--plumed=")==0){
      a.erase(0,a.find("=")+1);
      plumedFile=a;
      prefix="";
    } else if(a=="--plumed"){
      prefix="--plumed=";
    } else if(a.find("--timestep=")==0){
      a.erase(0,a.find("=")+1);
      double t;
      Tools::convert(a,t);
      timestep=real(t);
      prefix="";
    } else if(a=="--timestep"){
      prefix="--timestep=";
    } else if(a.find("--stride=")==0){
      a.erase(0,a.find("=")+1);
      Tools::convert(a,stride);
      prefix="";
    } else if(a=="--stride"){
      prefix="--stride=";
    } else if(a.find("--dump-forces=")==0){
      a.erase(0,a.find("=")+1);
      dumpforces=a;
      prefix="";
    } else if(a=="--dump-forces"){
      prefix="--dump-forces=";
    } else if(a.find("--dump-forces-fmt=")==0){
      a.erase(0,a.find("=")+1);
      dumpforcesFmt=a;
      prefix="";
    } else if(a=="--dump-forces-fmt"){
      prefix="--dump-forces-fmt=";
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
 "  [--help|-h]             : prints this help\n"
 "  [--help-debug]          : prints this help plus special options for debug\n"
 "  [--plumed FILE]         : plumed script file (default: plumed.dat)\n"
 "  [--timestep TS]         : timestep (default: 0.001) in picoseconds\n"
 "  [--stride ST]           : stride between frames (default: 1)\n"
 "  [--dump-forces FILE]    : dump forces on file FILE (default: do not dump)\n"
 "  [--dump-forces-fmt FMT] : dump forces on file FILE (default: %f)\n"
);
  if(printhelpdebug)
    fprintf(out,"%s",
 "Additional options for debug (only to be used in regtest):\n"
 "  [--debug-float]         : turns on the single precision version (to check float interface)\n"
);
    return 0;
  }


  if(trajectoryFile.length()==0){
    string msg="ERROR: please specify a trajectory";
    fprintf(stderr,"%s\n",msg.c_str());
    return 1;
  }

  Plumed p;
  int rr=sizeof(real);
  p.cmd("setRealPrecision",&rr);
  int checknatoms=0;
  int step=0;
  

  FILE* fp=fopen(trajectoryFile.c_str(),"r");

  FILE* fp_forces=NULL;
  if(dumpforces.length()>0){
    if(PlumedCommunicator::initialized() && pc.Get_size()>1){
      string n;
      Tools::convert(pc.Get_rank(),n);
      dumpforces+="."+n;
    }
    fp_forces=fopen(dumpforces.c_str(),"w");
  }
  
  std::string line;
  std::vector<real> coordinates;
  std::vector<real> forces;
  std::vector<real> masses;
  std::vector<real> cell;
  std::vector<real> virial;

  while(Tools::getline(fp,line)){

    int natoms;
    bool ok;
    sscanf(line.c_str(),"%d",&natoms);
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

    coordinates.assign(3*natoms,real(0.0));
    forces.assign(3*natoms,real(0.0));
    masses.assign(natoms,real(1.0));
    cell.assign(9,real(0.0));
    virial.assign(9,real(0.0));

    ok=Tools::getline(fp,line);
    plumed_massert(ok,"premature end of file");

    std::vector<std::string> words;
    words=Tools::getWords(line);
    std::vector<double> celld(9,0.0);
    if(words.size()==3){
      sscanf(line.c_str(),"%lf %lf %lf",&celld[0],&celld[4],&celld[8]);
    } else if(words.size()==9){
      sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
             &celld[0], &celld[1], &celld[2],
             &celld[3], &celld[4], &celld[5],
             &celld[6], &celld[7], &celld[8]);
    } else plumed_merror("needed box in second line");
    for(unsigned i=0;i<9;i++)cell[i]=real(celld[i]);
    for(int i=0;i<natoms;i++){
      ok=Tools::getline(fp,line);
      plumed_massert(ok,"premature end of file");
      char dummy[1000];
      double cc[3];
      std::sscanf(line.c_str(),"%s %lf %lf %lf",dummy,&cc[0],&cc[1],&cc[2]);
      coordinates[3*i]=real(cc[0]);
      coordinates[3*i+1]=real(cc[1]);
      coordinates[3*i+2]=real(cc[2]);
    }

   p.cmd("setForces",&forces[0]);
   p.cmd("setPositions",&coordinates[0]);
   p.cmd("setMasses",&masses[0]);
   p.cmd("setBox",&cell[0]);
   p.cmd("setVirial",&virial[0]);
   p.cmd("setStep",&step);
   p.cmd("calc");

// this is necessary as only processor zero is adding to the virial:
   pc.Bcast(&virial[0],9,0);

   if(fp_forces){
     fprintf(fp_forces,"%d\n",natoms);
// I use only a few digits for the forces since this is meant to be used
// with the regtests. Probably there should be an option for this
     string fmt=dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+"\n";
     fprintf(fp_forces,fmt.c_str(),virial[0],virial[4],virial[8]);
     fmt="X "+fmt;
     for(int i=0;i<natoms;i++)
       fprintf(fp_forces,fmt.c_str(),forces[3*i],forces[3*i+1],forces[3*i+2]);
   }

    step+=stride;
  }

  if(fp_forces) fclose(fp_forces);
  fclose(fp);
  
  return 0;
}

typedef CLToolDriver<double> Driver;
typedef CLToolDriver<float>  DriverSp;

PLUMED_REGISTER_CLTOOL(Driver,"driver")




}
