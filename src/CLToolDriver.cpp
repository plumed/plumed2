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
 by performing the actions described in the input file plumed.dat.  Actions that take the
stride keyword will be run for every frame in the trajectory.  The specific values of the 
STRIDE parameters in the input for PRINT,METAD,HISTOGRAM,etc will be ignored.
\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz
\endverbatim

The following command tells plumed to postprocess the trajectory contained in trajectory.xyz.
 by performing the actions described in the input file plumed.dat. Here though
--trajectory-stride is set equal to the frequency with which frames were output during the trajectory
and the --timestep is equal to the simulation timestep.  As such the STRIDE parameters in the plumed.dat
files are no longer ignored and any files output resemble those that would have been generated
had we run the calculation we are running with driver when the MD simulation was running.
\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --trajectory-stride 100 --timestep 0.001
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
  keys.add("compulsory","--timestep","1.0","the timestep that was used in the calculation that produced this trajectory in picoseconds");
  keys.add("compulsory","--stride","1","the frequency with which frames were output to this trajectory during the simulation");
  keys.add("atoms","--ixyz","the trajectory in xyz format");
  keys.add("optional","--length-units","units for length, either as a string or a number");
  keys.add("optional","--dump-forces","dump the forces on a file");
  keys.add("optional","--dump-forces-fmt","( default=%%f ) the format to use to dump the forces");
  keys.add("optional","--pdb","provides a pdb with masses and charges");
  keys.add("optional","--box","comma-separated box dimensions (3 for orthorombic, 9 for generic)");
  keys.add("hidden","--debug-float","turns on the single precision version (to check float interface)");
  keys.add("hidden","--debug-dd","use a fake domain decomposition");
  keys.add("hidden","--debug-pd","use a fake particle decomposition");
  keys.add("hidden","--debug-grex","use a fake gromacs-like replica exchange");
  keys.add("hidden","--debug-grex-log","log file for debug=grex");
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
  int multi=1;
  FILE*multi_log=NULL;
  bool debug_grex=parse("--debug-grex",fakein);
  PlumedCommunicator intracomm;
  PlumedCommunicator intercomm;
  if(debug_grex){
    Tools::convert(fakein,multi);
    int ntot=pc.Get_size();
    int nintra=ntot/multi;
    plumed_massert(multi*nintra==ntot,"xxx");
    pc.Split(pc.Get_rank()/nintra,pc.Get_rank(),intracomm);
    pc.Split(pc.Get_rank()%nintra,pc.Get_rank(),intercomm);
    string n; Tools::convert(intercomm.Get_rank(),n);
    string file;
    parse("--debug-grex-log",file);
    if(file.length()>0){
      file+="."+n;
      multi_log=fopen(file.c_str(),"w");
    }
  } else {
    intracomm.Set_comm(pc.Get_comm());
  }

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

  string pbc_cli_list; parse("--box",pbc_cli_list);
  bool pbc_cli_given=false;
  vector<double> pbc_cli_box(9,0.0);
  if(pbc_cli_list.length()>0) {
    pbc_cli_given=true;
    vector<string> words=Tools::getWords(pbc_cli_list,",");
    if(words.size()==3){
      for(int i=0;i<3;i++) sscanf(words[i].c_str(),"%lf",&(pbc_cli_box[4*i]));
    } else if(words.size()==9) {
      for(int i=0;i<9;i++) sscanf(words[i].c_str(),"%lf",&(pbc_cli_box[i]));
    } else {
      string msg="ERROR: cannot parse command-line box "+pbc_cli_list;
      fprintf(stderr,"%s\n",msg.c_str());
      return 1;
    }

  }
  
  plumed_massert(!(debug_dd&debug_pd),"cannot use debug-dd and debug-pd at the same time");
  if(debug_pd || debug_dd) plumed_massert(PlumedCommunicator::initialized(),"needs mpi for debug-pd");

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

  if(debug_grex){
    string n;
    Tools::convert(intercomm.Get_rank(),n);
    trajectoryFile+="."+n;
  }

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
// variables to test random decomposition (=domain decomposition)
  std::vector<int>  dd_gatindex;
  std::vector<int>  dd_g2l;
  std::vector<real> dd_masses;
  std::vector<real> dd_charges;
  std::vector<real> dd_forces;
  std::vector<real> dd_coordinates;
  int dd_nlocal;
// random stream to choose decompositions
  Random rnd;

  while(Tools::getline(fp,line)){

    int natoms;
    bool ok;
    bool first_step=false;
    sscanf(line.c_str(),"%d",&natoms);
    if(checknatoms==0){
      checknatoms=natoms;
      if(PlumedCommunicator::initialized()){
        if(multi>1){
          if(intracomm.Get_rank()==0) p.cmd("GREX setMPIIntercomm",&intercomm.Get_comm());
          p.cmd("GREX setMPIIntracomm",&intracomm.Get_comm());
          p.cmd("GREX init");
        }
        p.cmd("setMPIComm",&intracomm.Get_comm());
      }
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

    }
    plumed_massert(checknatoms==natoms,"number of atom changed");

    coordinates.assign(3*natoms,real(0.0));
    forces.assign(3*natoms,real(0.0));
    cell.assign(9,real(0.0));
    virial.assign(9,real(0.0));

    if( first_step || rnd.U01()>0.5){
      if(debug_pd){
        int npe=intracomm.Get_size();
        vector<int> loc(npe,0);
        vector<int> start(npe,0);
        for(int i=0;i<npe-1;i++){
          int cc=(natoms*2*rnd.U01())/npe;
          if(start[i]+cc>natoms) cc=natoms-start[i];
          loc[i]=cc;
          start[i+1]=start[i]+loc[i];
        }
        loc[npe-1]=natoms-start[npe-1];
        intracomm.Bcast(&loc[0],npe,0);
        intracomm.Bcast(&start[0],npe,0);
        pd_nlocal=loc[intracomm.Get_rank()];
        pd_start=start[intracomm.Get_rank()];
        if(intracomm.Get_rank()==0){
          fprintf(out,"\nDRIVER: Reassigning particle decomposition\n");
          fprintf(out,"DRIVER: "); for(int i=0;i<npe;i++) fprintf(out,"%d ",loc[i]); printf("\n");
          fprintf(out,"DRIVER: "); for(int i=0;i<npe;i++) fprintf(out,"%d ",start[i]); printf("\n");
        }
        p.cmd("setAtomsNlocal",&pd_nlocal);
        p.cmd("setAtomsContiguous",&pd_start);
      } else if(debug_dd){
        int npe=intracomm.Get_size();
        int rank=intracomm.Get_rank();
        dd_charges.assign(natoms,0.0);
        dd_masses.assign(natoms,0.0);
        dd_gatindex.assign(natoms,-1);
        dd_g2l.assign(natoms,-1);
        dd_coordinates.assign(3*natoms,0.0);
        dd_forces.assign(3*natoms,0.0);
        dd_nlocal=0;
        for(int i=0;i<natoms;++i){
          double r=rnd.U01()*npe;
          int n; for(n=0;n<npe;n++) if(n+1>r)break;
          plumed_assert(n<npe);
          if(n==rank){
            dd_gatindex[dd_nlocal]=i;
            dd_g2l[i]=dd_nlocal;
            dd_charges[dd_nlocal]=charges[i];
            dd_masses[dd_nlocal]=masses[i];
            dd_nlocal++;
          }
        }
        if(intracomm.Get_rank()==0){
          fprintf(out,"\nDRIVER: Reassigning particle decomposition\n");
        }
        p.cmd("setAtomsNlocal",&dd_nlocal);
        p.cmd("setAtomsGatindex",&dd_gatindex[0]);
      }
    }

    ok=Tools::getline(fp,line);
    plumed_massert(ok,"premature end of file");

    std::vector<double> celld(9,0.0);
    if(pbc_cli_given==false) {
      std::vector<std::string> words;
      words=Tools::getWords(line);
      if(words.size()==3){
	sscanf(line.c_str(),"%lf %lf %lf",&celld[0],&celld[4],&celld[8]);
      } else if(words.size()==9){
	sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &celld[0], &celld[1], &celld[2],
	       &celld[3], &celld[4], &celld[5],
	       &celld[6], &celld[7], &celld[8]);
      } else plumed_merror("needed box in second line");
    } else {			// from command line
      celld=pbc_cli_box;
    }
    for(unsigned i=0;i<9;i++)cell[i]=real(celld[i]);

    // Read coordinates
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
      if(debug_dd){
        for(int i=0;i<dd_nlocal;++i){
          int kk=dd_gatindex[i];
          dd_coordinates[3*i+0]=coordinates[3*kk+0];
          dd_coordinates[3*i+1]=coordinates[3*kk+1];
          dd_coordinates[3*i+2]=coordinates[3*kk+2];
        }
      }
    }

   if(debug_dd){
     p.cmd("setForces",&dd_forces[0]);
     p.cmd("setPositions",&dd_coordinates[0]);
     p.cmd("setMasses",&dd_masses[0]);
     p.cmd("setCharges",&dd_charges[0]);
   } else {
     p.cmd("setForces",&forces[3*pd_start]);
     p.cmd("setPositions",&coordinates[3*pd_start]);
     p.cmd("setMasses",&masses[pd_start]);
     p.cmd("setCharges",&charges[pd_start]);
   }
   p.cmd("setBox",&cell[0]);
   p.cmd("setVirial",&virial[0]);
   p.cmd("setStep",&step);
   p.cmd("calc");

// this is necessary as only processor zero is adding to the virial:
   intracomm.Bcast(&virial[0],9,0);
   if(debug_pd) intracomm.Sum(&forces[0],natoms*3);
   if(debug_dd){
     for(int i=0;i<dd_nlocal;i++){
       forces[3*dd_gatindex[i]+0]=dd_forces[3*i+0];
       forces[3*dd_gatindex[i]+1]=dd_forces[3*i+1];
       forces[3*dd_gatindex[i]+2]=dd_forces[3*i+2];
     }
     dd_forces.assign(3*natoms,0.0);
     intracomm.Sum(&forces[0],natoms*3);
   }
   int multi_stride=2;
   if(multi>1 &&step%multi_stride==0){
     p.cmd("GREX savePositions");
     if(intracomm.Get_rank()>0){
       p.cmd("GREX prepare");
     } else {
       int r=intercomm.Get_rank();
       int n=intercomm.Get_size();
       int partner=r+(2*((r+step/multi_stride)%2))-1;
       if(partner<0)partner=0;
       if(partner>=n) partner=n-1;
       p.cmd("GREX setPartner",&partner);
       p.cmd("GREX calculate");
       p.cmd("GREX shareAllDeltaBias");
       for(int i=0;i<n;i++){
         string s; Tools::convert(i,s);
         real a; s="GREX getDeltaBias "+s; p.cmd(s.c_str(),&a);
         if(multi_log) fprintf(multi_log," %f",a);
       }
       if(multi_log) fprintf(multi_log,"\n");
     }
   }


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
  if(multi_log) fclose(multi_log);
  
  return 0;
}

typedef CLToolDriver<double> CLToolDriverDouble;
typedef CLToolDriver<float> CLToolDriverFloat;
PLUMED_REGISTER_CLTOOL(CLToolDriverDouble,"driver")
PLUMED_REGISTER_CLTOOL(CLToolDriverFloat,"driver-float")




}
