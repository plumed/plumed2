/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "CLTool.h"
#include "CLToolRegister.h"
#include "Tools.h"
#include "Plumed.h"
#include "PlumedCommunicator.h"
#include "Random.h"
#include <cstdio>
#include <string>
#include <vector>
#include "Units.h"
#include "PDB.h"

using namespace std;

namespace PLMD {

//+PLUMEDOC TOOLS driver
/*
driver is a tool that allows one to to use plumed to post-process an existing trajectory.

The input to driver is specified using the command line arguments described below.

\par Examples

The following command tells plumed to postprocess the trajectory contained in trajectory.xyz
 by calcualting the actions described in the input file plumed.dat.
\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz
\endverbatim

*/
//+ENDPLUMEDOC

template<typename real>
class CLToolDriver : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  CLToolDriver(const CLToolOptions& co );
  int main(FILE* in,FILE*out,PlumedCommunicator& pc);
  string description()const;
};

template<typename real>
void CLToolDriver<real>::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys );
  keys.addFlag("--help-debug",false,"print special options that can be used to create regtests");
  keys.add("compulsory","--plumed","plumed.dat","specify the name of the plumed input file");
  keys.add("compulsory","--timestep","0.001","the timestep for the trajectory in picoseconds");
  keys.add("compulsory","--stride","1","stride between frames on which to do the calculation");
  keys.add("atoms","--ixyz","the trajectory in xyz format");
  keys.add("optional","--length-units","units for length, either as a string or a number");
  keys.add("optional","--dump-forces","dump the forces on a file");
  keys.add("optional","--dump-forces-fmt","( default=%%f ) the format to use to dump the forces");
  keys.add("optional","--pdb","provides a pdb with masses and charges");
  keys.add("hidden","--debug-float","turns on the single precision version (to check float interface)");
  keys.add("hidden","--debug-dd","use a fake domain decomposition");
  keys.add("hidden","--debug-pd","use a fake particle decomposition");
}

template<typename real>
CLToolDriver<real>::CLToolDriver(const CLToolOptions& co ):
CLTool(co)
{
 inputdata=commandline;
}

template<typename real>
string CLToolDriver<real>::description()const{ return "analyze trajectories with plumed"; }

template<>
string CLToolDriver<float>::description()const{ return "analyze trajectories with plumed (single precision version)"; }


template<typename real>
int CLToolDriver<real>::main(FILE* in,FILE*out,PlumedCommunicator& pc){

  Units units;
  PDB pdb;

// Parse everything
  bool printhelpdebug; parseFlag("--help-debug",printhelpdebug);
  if( printhelpdebug ){
      fprintf(out,"%s",
         "Additional options for debug (only to be used in regtest):\n"
         "  [--debug-float]         : turns on the single precision version (to check float interface)\n"
         "  [--debug-dd]            : use a fake domain decomposition\n"
         "  [--debug-pd]            : use a fake particle decomposition\n"
      );
      return 0;
  }
  std::string fakein; 
  bool debugfloat=parse("--debug-float",fakein);
  if(debugfloat && sizeof(real)!=sizeof(float)){
      CLTool* cl=cltoolRegister().create(CLToolOptions("driver-float"));    //new CLToolDriver<float>(*this);
      cl->inputData=this->inputData; 
      int ret=cl->main(in,out,pc);
      delete cl;
      return ret;
  }

  bool debug_pd=parse("--debug-pd",fakein);
  bool debug_dd=parse("--debug-dd",fakein);
// Read the plumed input file name  
  string plumedFile; parse("--plumed",plumedFile);
// the timestep
  double t; parse("--timestep",t);
  real timestep=real(t);
// the stride
  unsigned stride; parse("--stride",stride);
// are we writing forces
  string dumpforces(""), dumpforcesFmt("%f");; 
  parse("--dump-forces",dumpforces);
  if(dumpforces!="") parse("--dump-forces-fmt",dumpforcesFmt);

// Read in an xyz file
  string trajectoryFile("");
  std::string traj_xyz; parse("--ixyz",traj_xyz);
  if(traj_xyz.length()>0 && trajectoryFile.length()==0) trajectoryFile=traj_xyz;
  if(trajectoryFile.length()==0){
    fprintf(out,"ERROR: missing trajectory data\n"); 
    return 0;
  }
  string lengthUnits(""); parse("--length-units",lengthUnits);
  if(lengthUnits.length()>0) units.setLength(lengthUnits);

  string pdbfile; parse("--pdb",pdbfile);
  if(pdbfile.length()>0){
    bool check=pdb.read(pdbfile,false,1.0);
    plumed_massert(check,"error reading pdb file");
  }

  if(debug_dd) plumed_merror("debug_dd not yet implemented");
  plumed_massert(!(debug_dd&debug_pd),"cannot use debug-dd and debug-pd at the same time");
  if(debug_pd) plumed_massert(PlumedCommunicator::initialized(),"needs mpi for debug-pd");

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

  FILE* fp;
  if (trajectoryFile=="-") 
    fp=in;
  else {
    fp=fopen(trajectoryFile.c_str(),"r");
    if(!fp){
      string msg="ERROR: Error opening XYZ file "+trajectoryFile;
      fprintf(stderr,"%s\n",msg.c_str());
      return 1;
    }
  }
    

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
  std::vector<real> charges;
  std::vector<real> cell;
  std::vector<real> virial;

// variables to test particle decomposition
  int pd_nlocal;
  int pd_start;
// random stream to choose decompositions
  Random rnd;

  while(Tools::getline(fp,line)){

    int natoms;
    bool ok;
    bool first_step=false;
    sscanf(line.c_str(),"%d",&natoms);
    if(checknatoms==0){
      checknatoms=natoms;
      if(PlumedCommunicator::initialized()) p.cmd("setMPIComm",&pc.Get_comm());
      p.cmd("setMDLengthUnits",&units.getLength());
      p.cmd("setNatoms",&natoms);
      p.cmd("setMDEngine","driver");
      p.cmd("setTimestep",&timestep);
      p.cmd("setPlumedDat",plumedFile.c_str());
      p.cmd("setLog",out);
      p.cmd("init");
      pd_nlocal=natoms;
      pd_start=0;
      first_step=true;
    }
    plumed_massert(checknatoms==natoms,"number of atom changed");

    coordinates.assign(3*natoms,real(0.0));
    forces.assign(3*natoms,real(0.0));
    masses.assign(natoms,real(1.0));
    charges.assign(natoms,real(0.0));
    if(pdbfile.length()>0){
      for(unsigned i=0;i<pdb.size();++i){
        AtomNumber an=pdb.getAtomNumbers()[i];
        unsigned index=an.index();
        plumed_massert(index<unsigned(natoms),"atom index in pdb exceeds the number of atoms in trajectory");
        masses[index]=pdb.getOccupancy()[i];
        charges[index]=pdb.getBeta()[i];
      }
    }
    cell.assign(9,real(0.0));
    virial.assign(9,real(0.0));

    if(debug_pd && ( first_step || rnd.U01()>0.5)){
      int npe=pc.Get_size();
      vector<int> loc(npe,0);
      vector<int> start(npe,0);
      for(int i=0;i<npe-1;i++){
        int cc=(natoms*2*rnd.U01())/npe;
        if(start[i]+cc>natoms) cc=natoms-start[i];
        loc[i]=cc;
        start[i+1]=start[i]+loc[i];
      }
      loc[npe-1]=natoms-start[npe-1];
      pc.Bcast(&loc[0],npe,0);
      pc.Bcast(&start[0],npe,0);
      pd_nlocal=loc[pc.Get_rank()];
      pd_start=start[pc.Get_rank()];
      if(pc.Get_rank()==0){
        fprintf(out,"\nDRIVER: Reassigning particle decomposition\n");
        fprintf(out,"DRIVER: "); for(int i=0;i<npe;i++) fprintf(out,"%d ",loc[i]); printf("\n");
        fprintf(out,"DRIVER: "); for(int i=0;i<npe;i++) fprintf(out,"%d ",start[i]); printf("\n");
      }
      p.cmd("setAtomsNlocal",&pd_nlocal);
      p.cmd("setAtomsContiguous",&pd_start);
    }

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
      if(!debug_pd || ( i>=pd_start && i<pd_start+pd_nlocal) ){
        coordinates[3*i]=real(cc[0]);
        coordinates[3*i+1]=real(cc[1]);
        coordinates[3*i+2]=real(cc[2]);
      }
    }

   p.cmd("setForces",&forces[3*pd_start]);
   p.cmd("setPositions",&coordinates[3*pd_start]);
   p.cmd("setMasses",&masses[3*pd_start]);
   p.cmd("setCharges",&charges[3*pd_start]);
   p.cmd("setBox",&cell[0]);
   p.cmd("setVirial",&virial[0]);
   p.cmd("setStep",&step);
   p.cmd("calc");

// this is necessary as only processor zero is adding to the virial:
   pc.Bcast(&virial[0],9,0);
   if(debug_pd) pc.Sum(&forces[0],natoms*3);

   if(fp_forces){
     fprintf(fp_forces,"%d\n",natoms);
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

typedef CLToolDriver<double> CLToolDriverDouble;
typedef CLToolDriver<float> CLToolDriverFloat;
PLUMED_REGISTER_CLTOOL(CLToolDriverDouble,"driver")
PLUMED_REGISTER_CLTOOL(CLToolDriverFloat,"driver-float")




}
