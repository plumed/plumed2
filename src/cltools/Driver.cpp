/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

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
#include "tools/Tools.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/Random.h"
#include "tools/Pbc.h"
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include <memory>
#include "tools/Units.h"
#include "tools/PDB.h"
#include "tools/FileBase.h"
#include "tools/IFile.h"

// when using molfile plugin
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
#ifndef __PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS
/* Use the internal ones. Alternatively:
 *    ifeq (,$(findstring __PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS,$(CPPFLAGS)))
 *    CPPFLAGS+=-I../molfile
 */
#include "molfile/libmolfile_plugin.h"
#include "molfile/molfile_plugin.h"
using namespace PLMD::molfile;
#else
#include <libmolfile_plugin.h>
#include <molfile_plugin.h>
#endif
#endif

#ifdef __PLUMED_HAS_XDRFILE
#include <xdrfile/xdrfile_trr.h>
#include <xdrfile/xdrfile_xtc.h>
#endif

using namespace std;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS driver-float
/*
Equivalent to \ref driver, but using single precision reals.

The purpose of this tool is just to test what PLUMED does when linked from
a single precision code.

\par Examples

\verbatim
plumed driver-float --plumed plumed.dat --ixyz trajectory.xyz
\endverbatim

See also examples in \ref driver

*/
//+ENDPLUMEDOC
//


//+PLUMEDOC TOOLS driver
/*
driver is a tool that allows one to to use plumed to post-process an existing trajectory.

The input to driver is specified using the command line arguments described below.

In addition, you can use the special \subpage READ command inside your plumed input
to read in colvar files that were generated during your MD simulation.  The values
read in can then be treated like calculated colvars.

\warning
Notice that by default the driver has no knowledge about the masses and charges
of your atoms! Thus, if you want to compute quantities depending charges (e.g. \ref DHENERGY)
or masses (e.g. \ref COM) you should pass the proper information to the driver.
You can do it either with the --pdb option or with the --mc option. The latter
will read a file produced by \ref DUMPMASSCHARGE .


\par Examples

The following command tells plumed to post process the trajectory contained in `trajectory.xyz`
 by performing the actions described in the input file `plumed.dat`.  If an action that takes the
stride keyword is given a stride equal to \f$n\f$ then it will be performed only on every \f$n\f$th
frames in the trajectory file.
\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz
\endverbatim

Notice that `xyz` files are expected to be in internal PLUMED units, that is by default nm.
You can change this behavior by using the `--length-units` option:
\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --length-units A
\endverbatim
The strings accepted by the `--length-units` options are the same ones accepted by the \ref UNITS action.
Other file formats typically have their default coordinates (e.g., `gro` files are always in nm)
and it thus should not be necessary to use the `--length-units` option. Additionally,
consider that the units used by the `driver` might be different by the units used in the PLUMED input
file `plumed.dat`. For instance consider the command:
\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --length-units A
\endverbatim
where `plumed.dat` is
\plumedfile
# no explicit UNITS action here
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=colvar
\endplumedfile
In this case, the driver reads the `xyz` file assuming it to contain coordinates in Angstrom units.
However, the resulting `colvar` file contains a distance expressed in nm.

The following command tells plumed to post process the trajectory contained in trajectory.xyz.
 by performing the actions described in the input file plumed.dat.
\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --trajectory-stride 100 --timestep 0.001
\endverbatim
Here though
`--trajectory-stride` is set equal to the frequency with which frames were output during the trajectory
and the `--timestep` is equal to the simulation timestep.  As such the `STRIDE` parameters in the `plumed.dat`
files are referred to the original timestep and any files output resemble those that would have been generated
had we run the calculation we are running with driver when the MD simulation was running.

PLUMED can read xyz files (in PLUMED units) and gro files (in nm). In addition,
PLUMED includes by default support for a
subset of the trajectory file formats supported by VMD, e.g. xtc and dcd:

\verbatim
plumed driver --plumed plumed.dat --pdb diala.pdb --mf_xtc traj.xtc --trajectory-stride 100 --timestep 0.001
\endverbatim

where `--mf_` prefixes the extension of one of the accepted molfile plugin format.
If PLUMED has been \ref Installation "installed" with full molfile support, other formats will be available.
Just type `plumed driver --help` to see which plugins are available.

Molfile plugin require periodic cell to be triangular (i.e. first vector oriented along x and
second vector in xy plane). This is true for many MD codes. However, it could be false
if you rotate the coordinates in your trajectory before reading them in the driver.
Also notice that some formats (e.g. amber crd) do not specify atom number. In this case you can use
the `--natoms` option:
\verbatim
plumed driver --plumed plumed.dat --imf_crd trajectory.crd --natoms 128
\endverbatim

Check the available molfile plugins and limitations at [this link](http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/).

Additionally, you can use the xdrfile implementation of xtc and trr. To this aim, just
download and install properly the xdrfile library (see [this link](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library)).
If the xdrfile library is installed properly the PLUMED configure script should be able to
detect it and enable it.
Notice that the xdrfile implementation of xtc and trr
is more robust than the molfile one, since it provides support for generic cell shapes.
In addition, it allows \ref DUMPATOMS to write compressed xtc files.


*/
//+ENDPLUMEDOC
//

#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
static vector<molfile_plugin_t *> plugins;
static map <string, unsigned> pluginmap;
static int register_cb(void *v, vmdplugin_t *p) {
  //const char *key = p->name;
  const auto ret = pluginmap.insert ( std::pair<string,unsigned>(string(p->name),plugins.size()) );
  if (ret.second==false) {
    //cerr<<"MOLFILE: found duplicate plugin for "<<key<<" : not inserted "<<endl;
  } else {
    //cerr<<"MOLFILE: loading plugin "<<key<<" number "<<plugins.size()-1<<endl;
    plugins.push_back(reinterpret_cast<molfile_plugin_t *>(p));
  }
  return VMDPLUGIN_SUCCESS;
}
#endif

template<typename real>
class Driver : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit Driver(const CLToolOptions& co );
  int main(FILE* in,FILE*out,Communicator& pc) override;
  void evaluateNumericalDerivatives( const long int& step, PlumedMain& p, const std::vector<real>& coordinates,
                                     const std::vector<real>& masses, const std::vector<real>& charges,
                                     std::vector<real>& cell, const double& base, std::vector<real>& numder );
  string description()const override;
};

template<typename real>
void Driver<real>::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys ); keys.isDriver();
  keys.addFlag("--help-debug",false,"print special options that can be used to create regtests");
  keys.add("compulsory","--plumed","plumed.dat","specify the name of the plumed input file");
  keys.add("compulsory","--timestep","1.0","the timestep that was used in the calculation that produced this trajectory in picoseconds");
  keys.add("compulsory","--trajectory-stride","1","the frequency with which frames were output to this trajectory during the simulation"
#ifdef __PLUMED_HAS_XDRFILE
           " (0 means that the number of the step is read from the trajectory file,"
           " currently working only for xtc/trr files read with --ixtc/--trr)"
#endif
          );
  keys.add("compulsory","--multi","0","set number of replicas for multi environment (needs MPI)");
  keys.addFlag("--noatoms",false,"don't read in a trajectory.  Just use colvar files as specified in plumed.dat");
  keys.addFlag("--parse-only",false,"read the plumed input file and stop");
  keys.add("atoms","--ixyz","the trajectory in xyz format");
  keys.add("atoms","--igro","the trajectory in gro format");
  keys.add("atoms","--idlp4","the trajectory in DL_POLY_4 format");
#ifdef __PLUMED_HAS_XDRFILE
  keys.add("atoms","--ixtc","the trajectory in xtc format (xdrfile implementation)");
  keys.add("atoms","--itrr","the trajectory in trr format (xdrfile implementation)");
#endif
  keys.add("optional","--length-units","units for length, either as a string or a number");
  keys.add("optional","--mass-units","units for mass in pdb and mc file, either as a string or a number");
  keys.add("optional","--charge-units","units for charge in pdb and mc file, either as a string or a number");
  keys.add("optional","--kt","set \\f$k_B T\\f$, it will not be necessary to specify temperature in input file");
  keys.add("optional","--dump-forces","dump the forces on a file");
  keys.add("optional","--dump-forces-fmt","( default=%%f ) the format to use to dump the forces");
  keys.addFlag("--dump-full-virial",false,"with --dump-forces, it dumps the 9 components of the virial");
  keys.add("optional","--pdb","provides a pdb with masses and charges");
  keys.add("optional","--mc","provides a file with masses and charges as produced with DUMPMASSCHARGE");
  keys.add("optional","--box","comma-separated box dimensions (3 for orthorhombic, 9 for generic)");
  keys.add("optional","--natoms","provides number of atoms - only used if file format does not contain number of atoms");
  keys.add("optional","--initial-step","provides a number for the initial step, default is 0");
  keys.add("optional","--debug-forces","output a file containing the forces due to the bias evaluated using numerical derivatives "
           "and using the analytical derivatives implemented in plumed");
  keys.add("hidden","--debug-float","[yes/no] turns on the single precision version (to check float interface)");
  keys.add("hidden","--debug-dd","[yes/no] use a fake domain decomposition");
  keys.add("hidden","--debug-pd","[yes/no] use a fake particle decomposition");
  keys.add("hidden","--debug-grex","use a fake gromacs-like replica exchange, specify exchange stride");
  keys.add("hidden","--debug-grex-log","log file for debug=grex");
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
  MOLFILE_INIT_ALL
  MOLFILE_REGISTER_ALL(NULL, register_cb)
  for(unsigned i=0; i<plugins.size(); i++) {
    string kk="--mf_"+string(plugins[i]->name);
    string mm=" molfile: the trajectory in "+string(plugins[i]->name)+" format " ;
    keys.add("atoms",kk,mm);
  }
#endif
}
template<typename real>
Driver<real>::Driver(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}
template<typename real>
string Driver<real>::description()const { return "analyze trajectories with plumed"; }

template<typename real>
int Driver<real>::main(FILE* in,FILE*out,Communicator& pc) {

  Units units;
  PDB pdb;

// Parse everything
  bool printhelpdebug; parseFlag("--help-debug",printhelpdebug);
  if( printhelpdebug ) {
    fprintf(out,"%s",
            "Additional options for debug (only to be used in regtest):\n"
            "  [--debug-float yes]     : turns on the single precision version (to check float interface)\n"
            "  [--debug-dd yes]        : use a fake domain decomposition\n"
            "  [--debug-pd yes]        : use a fake particle decomposition\n"
           );
    return 0;
  }
  // Are we reading trajectory data
  bool noatoms; parseFlag("--noatoms",noatoms);
  bool parseOnly; parseFlag("--parse-only",parseOnly);

  std::string fakein;
  bool debug_float=false;
  fakein="";
  if(parse("--debug-float",fakein)) {
    if(fakein=="yes") debug_float=true;
    else if(fakein=="no") debug_float=false;
    else error("--debug-float should have argument yes or no");
  }

  if(debug_float && sizeof(real)!=sizeof(float)) {
    auto cl=cltoolRegister().create(CLToolOptions("driver-float"));
    cl->setInputData(this->getInputData());
    int ret=cl->main(in,out,pc);
    return ret;
  }

  bool debug_pd=false;
  fakein="";
  if(parse("--debug-pd",fakein)) {
    if(fakein=="yes") debug_pd=true;
    else if(fakein=="no") debug_pd=false;
    else error("--debug-pd should have argument yes or no");
  }
  if(debug_pd) fprintf(out,"DEBUGGING PARTICLE DECOMPOSITION\n");

  bool debug_dd=false;
  fakein="";
  if(parse("--debug-dd",fakein)) {
    if(fakein=="yes") debug_dd=true;
    else if(fakein=="no") debug_dd=false;
    else error("--debug-dd should have argument yes or no");
  }
  if(debug_dd) fprintf(out,"DEBUGGING DOMAIN DECOMPOSITION\n");

  if( debug_pd || debug_dd ) {
    if(noatoms) error("cannot debug without atoms");
  }

// set up for multi replica driver:
  int multi=0;
  parse("--multi",multi);
  Communicator intracomm;
  Communicator intercomm;
  if(multi) {
    int ntot=pc.Get_size();
    int nintra=ntot/multi;
    if(multi*nintra!=ntot) error("invalid number of processes for multi environment");
    pc.Split(pc.Get_rank()/nintra,pc.Get_rank(),intracomm);
    pc.Split(pc.Get_rank()%nintra,pc.Get_rank(),intercomm);
  } else {
    intracomm.Set_comm(pc.Get_comm());
  }

// set up for debug replica exchange:
  bool debug_grex=parse("--debug-grex",fakein);
  int  grex_stride=0;
  FILE*grex_log=NULL;
  if(debug_grex) {
    if(noatoms) error("must have atoms to debug_grex");
    if(multi<2)  error("--debug_grex needs --multi with at least two replicas");
    Tools::convert(fakein,grex_stride);
    string n; Tools::convert(intercomm.Get_rank(),n);
    string file;
    parse("--debug-grex-log",file);
    if(file.length()>0) {
      file+="."+n;
      grex_log=fopen(file.c_str(),"w");
    }
  }

// Read the plumed input file name
  string plumedFile; parse("--plumed",plumedFile);
// the timestep
  double t; parse("--timestep",t);
  real timestep=real(t);
// the stride
  unsigned stride; parse("--trajectory-stride",stride);
// are we writing forces
  string dumpforces(""), debugforces(""), dumpforcesFmt("%f");;
  bool dumpfullvirial=false;
  if(!noatoms) {
    parse("--dump-forces",dumpforces);
    parse("--debug-forces",debugforces);
  }
  if(dumpforces!="" || debugforces!="" ) parse("--dump-forces-fmt",dumpforcesFmt);
  if(dumpforces!="") parseFlag("--dump-full-virial",dumpfullvirial);
  if( debugforces!="" && (debug_dd || debug_pd) ) error("cannot debug forces and domain/particle decomposition at same time");
  if( debugforces!="" && sizeof(real)!=sizeof(double) ) error("cannot debug forces in single precision mode");

  real kt=-1.0;
  parse("--kt",kt);
  string trajectory_fmt;

  bool use_molfile=false;
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
  molfile_plugin_t *api=NULL;
  void *h_in=NULL;
  molfile_timestep_t ts_in; // this is the structure that has the timestep
// a std::unique_ptr<float> with the same scope as ts_in
// it is necessary in order to store the pointer to ts_in.coords
  std::unique_ptr<float[]> ts_in_coords;
  ts_in.coords=ts_in_coords.get();
  ts_in.A=-1; // we use this to check whether cell is provided or not
#endif

// Read in an xyz file
  string trajectoryFile(""), pdbfile(""), mcfile("");
  bool pbc_cli_given=false; vector<double> pbc_cli_box(9,0.0);
  int command_line_natoms=-1;

  if(!noatoms) {
    std::string traj_xyz; parse("--ixyz",traj_xyz);
    std::string traj_gro; parse("--igro",traj_gro);
    std::string traj_dlp4; parse("--idlp4",traj_dlp4);
    std::string traj_xtc;
    std::string traj_trr;
#ifdef __PLUMED_HAS_XDRFILE
    parse("--ixtc",traj_xtc);
    parse("--itrr",traj_trr);
#endif
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
    for(unsigned i=0; i<plugins.size(); i++) {
      string molfile_key="--mf_"+string(plugins[i]->name);
      string traj_molfile;
      parse(molfile_key,traj_molfile);
      if(traj_molfile.length()>0) {
        fprintf(out,"\nDRIVER: Found molfile format trajectory %s with name %s\n",plugins[i]->name,traj_molfile.c_str());
        trajectoryFile=traj_molfile;
        trajectory_fmt=string(plugins[i]->name);
        use_molfile=true;
        api = plugins[i];
      }
    }
#endif
    { // check that only one fmt is specified
      int nn=0;
      if(traj_xyz.length()>0) nn++;
      if(traj_gro.length()>0) nn++;
      if(traj_dlp4.length()>0) nn++;
      if(traj_xtc.length()>0) nn++;
      if(traj_trr.length()>0) nn++;
      if(nn>1) {
        fprintf(stderr,"ERROR: cannot provide more than one trajectory file\n");
        if(grex_log)fclose(grex_log);
        return 1;
      }
    }
    if(traj_xyz.length()>0 && trajectoryFile.length()==0) {
      trajectoryFile=traj_xyz;
      trajectory_fmt="xyz";
    }
    if(traj_gro.length()>0 && trajectoryFile.length()==0) {
      trajectoryFile=traj_gro;
      trajectory_fmt="gro";
    }
    if(traj_dlp4.length()>0 && trajectoryFile.length()==0) {
      trajectoryFile=traj_dlp4;
      trajectory_fmt="dlp4";
    }
    if(traj_xtc.length()>0 && trajectoryFile.length()==0) {
      trajectoryFile=traj_xtc;
      trajectory_fmt="xdr-xtc";
    }
    if(traj_trr.length()>0 && trajectoryFile.length()==0) {
      trajectoryFile=traj_trr;
      trajectory_fmt="xdr-trr";
    }
    if(trajectoryFile.length()==0&&!parseOnly) {
      fprintf(stderr,"ERROR: missing trajectory data\n");
      if(grex_log)fclose(grex_log);
      return 1;
    }
    string lengthUnits(""); parse("--length-units",lengthUnits);
    if(lengthUnits.length()>0) units.setLength(lengthUnits);
    string chargeUnits(""); parse("--charge-units",chargeUnits);
    if(chargeUnits.length()>0) units.setCharge(chargeUnits);
    string massUnits(""); parse("--mass-units",massUnits);
    if(massUnits.length()>0) units.setMass(massUnits);

    parse("--pdb",pdbfile);
    if(pdbfile.length()>0) {
      bool check=pdb.read(pdbfile,false,1.0);
      if(!check) error("error reading pdb file");
    }

    parse("--mc",mcfile);

    string pbc_cli_list; parse("--box",pbc_cli_list);
    if(pbc_cli_list.length()>0) {
      pbc_cli_given=true;
      vector<string> words=Tools::getWords(pbc_cli_list,",");
      if(words.size()==3) {
        for(int i=0; i<3; i++) sscanf(words[i].c_str(),"%100lf",&(pbc_cli_box[4*i]));
      } else if(words.size()==9) {
        for(int i=0; i<9; i++) sscanf(words[i].c_str(),"%100lf",&(pbc_cli_box[i]));
      } else {
        string msg="ERROR: cannot parse command-line box "+pbc_cli_list;
        fprintf(stderr,"%s\n",msg.c_str());
        return 1;
      }

    }

    parse("--natoms",command_line_natoms);

  }


  if(debug_dd && debug_pd) error("cannot use debug-dd and debug-pd at the same time");
  if(debug_pd || debug_dd) {
    if( !Communicator::initialized() ) error("needs mpi for debug-pd");
  }

  PlumedMain p;
  int rr=sizeof(real);
  p.cmd("setRealPrecision",&rr);
  int checknatoms=-1;
  long int step=0;
  parse("--initial-step",step);

  if(Communicator::initialized()) {
    if(multi) {
      if(intracomm.Get_rank()==0) p.cmd("GREX setMPIIntercomm",&intercomm.Get_comm());
      p.cmd("GREX setMPIIntracomm",&intracomm.Get_comm());
      p.cmd("GREX init");
    }
    p.cmd("setMPIComm",&intracomm.Get_comm());
  }
  p.cmd("setMDLengthUnits",&units.getLength());
  p.cmd("setMDChargeUnits",&units.getCharge());
  p.cmd("setMDMassUnits",&units.getMass());
  p.cmd("setMDEngine","driver");
  p.cmd("setTimestep",&timestep);
  p.cmd("setPlumedDat",plumedFile.c_str());
  p.cmd("setLog",out);

  int natoms;
  int lvl=0;
  int pb=1;

  if(parseOnly) {
    if(command_line_natoms<0) error("--parseOnly requires setting the number of atoms with --natoms");
    natoms=command_line_natoms;
  }


  FILE* fp=NULL; FILE* fp_forces=NULL; OFile fp_dforces;
#ifdef __PLUMED_HAS_XDRFILE
  XDRFILE* xd=NULL;
#endif
  if(!noatoms&&!parseOnly) {
    if (trajectoryFile=="-")
      fp=in;
    else {
      if(multi) {
        string n;
        Tools::convert(intercomm.Get_rank(),n);
        std::string testfile=FileBase::appendSuffix(trajectoryFile,"."+n);
        FILE* tmp_fp=fopen(testfile.c_str(),"r");
        if(tmp_fp) { fclose(tmp_fp); trajectoryFile=testfile.c_str();}
      }
      if(use_molfile==true) {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
        h_in = api->open_file_read(trajectoryFile.c_str(), trajectory_fmt.c_str(), &natoms);
        if(natoms==MOLFILE_NUMATOMS_UNKNOWN) {
          if(command_line_natoms>=0) natoms=command_line_natoms;
          else error("this file format does not provide number of atoms; use --natoms on the command line");
        }
        ts_in_coords.reset(new float [3*natoms]);
        ts_in.coords = ts_in_coords.get();
#endif
      } else if(trajectory_fmt=="xdr-xtc" || trajectory_fmt=="xdr-trr") {
#ifdef __PLUMED_HAS_XDRFILE
        xd=xdrfile_open(trajectoryFile.c_str(),"r");
        if(!xd) {
          string msg="ERROR: Error opening trajectory file "+trajectoryFile;
          fprintf(stderr,"%s\n",msg.c_str());
          return 1;
        }
        if(trajectory_fmt=="xdr-xtc") read_xtc_natoms(&trajectoryFile[0],&natoms);
        if(trajectory_fmt=="xdr-trr") read_trr_natoms(&trajectoryFile[0],&natoms);
#endif
      } else {
        fp=fopen(trajectoryFile.c_str(),"r");
        if(!fp) {
          string msg="ERROR: Error opening trajectory file "+trajectoryFile;
          fprintf(stderr,"%s\n",msg.c_str());
          return 1;
        }
      }
    }
    if(dumpforces.length()>0) {
      if(Communicator::initialized() && pc.Get_size()>1) {
        string n;
        Tools::convert(pc.Get_rank(),n);
        dumpforces+="."+n;
      }
      fp_forces=fopen(dumpforces.c_str(),"w");
    }
    if(debugforces.length()>0) {
      if(Communicator::initialized() && pc.Get_size()>1) {
        string n;
        Tools::convert(pc.Get_rank(),n);
        debugforces+="."+n;
      }
      fp_dforces.open(debugforces);
    }
  }

  std::string line;
  std::vector<real> coordinates;
  std::vector<real> forces;
  std::vector<real> masses;
  std::vector<real> charges;
  std::vector<real> cell;
  std::vector<real> virial;
  std::vector<real> numder;

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

  if(trajectory_fmt=="dlp4") {
    if(!Tools::getline(fp,line)) error("error reading title");
    if(!Tools::getline(fp,line)) error("error reading atoms");
    sscanf(line.c_str(),"%d %d %d",&lvl,&pb,&natoms);

  }
  bool lstep=true;
  while(true) {
    if(!noatoms&&!parseOnly) {
      if(use_molfile==true) {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
        int rc;
        rc = api->read_next_timestep(h_in, natoms, &ts_in);
        if(rc==MOLFILE_EOF) {
          break;
        }
#endif
      } else if(trajectory_fmt=="xyz" || trajectory_fmt=="gro" || trajectory_fmt=="dlp4") {
        if(!Tools::getline(fp,line)) break;
      }
    }
    bool first_step=false;
    if(!noatoms&&!parseOnly) {
      if(use_molfile==false && (trajectory_fmt=="xyz" || trajectory_fmt=="gro")) {
        if(trajectory_fmt=="gro") if(!Tools::getline(fp,line)) error("premature end of trajectory file");
        sscanf(line.c_str(),"%100d",&natoms);
      }
      if(use_molfile==false && trajectory_fmt=="dlp4") {
        char xa[9];
        int xb,xc,xd;
        double t;
        sscanf(line.c_str(),"%8s %ld %d %d %d %lf",xa,&step,&xb,&xc,&xd,&t);
        timestep = real(t);
        if (lstep) {
          p.cmd("setTimestep",&timestep);
          lstep = false;
        }
      }
    }
    if(checknatoms<0 && !noatoms) {
      pd_nlocal=natoms;
      pd_start=0;
      first_step=true;
      masses.assign(natoms,std::numeric_limits<real>::quiet_NaN());
      charges.assign(natoms,std::numeric_limits<real>::quiet_NaN());
//case pdb: structure
      if(pdbfile.length()>0) {
        for(unsigned i=0; i<pdb.size(); ++i) {
          AtomNumber an=pdb.getAtomNumbers()[i];
          unsigned index=an.index();
          if( index>=unsigned(natoms) ) error("atom index in pdb exceeds the number of atoms in trajectory");
          masses[index]=pdb.getOccupancy()[i];
          charges[index]=pdb.getBeta()[i];
        }
      }
      if(mcfile.length()>0) {
        IFile ifile;
        ifile.open(mcfile);
        int index; double mass; double charge;
        while(ifile.scanField("index",index).scanField("mass",mass).scanField("charge",charge).scanField()) {
          masses[index]=mass;
          charges[index]=charge;
        }
      }
    } else if( checknatoms<0 && noatoms ) {
      natoms=0;
    }
    if( checknatoms<0 ) {
      if(kt>=0) {
        p.cmd("setKbT",&kt);
      }
      checknatoms=natoms;
      p.cmd("setNatoms",&natoms);
      p.cmd("init");
      if(parseOnly) break;
    }
    if(checknatoms!=natoms) {
      std::string stepstr; Tools::convert(step,stepstr);
      error("number of atoms in frame " + stepstr + " does not match number of atoms in first frame");
    }

    coordinates.assign(3*natoms,real(0.0));
    forces.assign(3*natoms,real(0.0));
    cell.assign(9,real(0.0));
    virial.assign(9,real(0.0));

    if( first_step || rnd.U01()>0.5) {
      if(debug_pd) {
        int npe=intracomm.Get_size();
        vector<int> loc(npe,0);
        vector<int> start(npe,0);
        for(int i=0; i<npe-1; i++) {
          int cc=(natoms*2*rnd.U01())/npe;
          if(start[i]+cc>natoms) cc=natoms-start[i];
          loc[i]=cc;
          start[i+1]=start[i]+loc[i];
        }
        loc[npe-1]=natoms-start[npe-1];
        intracomm.Bcast(loc,0);
        intracomm.Bcast(start,0);
        pd_nlocal=loc[intracomm.Get_rank()];
        pd_start=start[intracomm.Get_rank()];
        if(intracomm.Get_rank()==0) {
          fprintf(out,"\nDRIVER: Reassigning particle decomposition\n");
          fprintf(out,"DRIVER: "); for(int i=0; i<npe; i++) fprintf(out,"%d ",loc[i]); printf("\n");
          fprintf(out,"DRIVER: "); for(int i=0; i<npe; i++) fprintf(out,"%d ",start[i]); printf("\n");
        }
        p.cmd("setAtomsNlocal",&pd_nlocal);
        p.cmd("setAtomsContiguous",&pd_start);
      } else if(debug_dd) {
        int npe=intracomm.Get_size();
        int rank=intracomm.Get_rank();
        dd_charges.assign(natoms,0.0);
        dd_masses.assign(natoms,0.0);
        dd_gatindex.assign(natoms,-1);
        dd_g2l.assign(natoms,-1);
        dd_coordinates.assign(3*natoms,0.0);
        dd_forces.assign(3*natoms,0.0);
        dd_nlocal=0;
        for(int i=0; i<natoms; ++i) {
          double r=rnd.U01()*npe;
          int n; for(n=0; n<npe; n++) if(n+1>r)break;
          plumed_assert(n<npe);
          if(n==rank) {
            dd_gatindex[dd_nlocal]=i;
            dd_g2l[i]=dd_nlocal;
            dd_charges[dd_nlocal]=charges[i];
            dd_masses[dd_nlocal]=masses[i];
            dd_nlocal++;
          }
        }
        if(intracomm.Get_rank()==0) {
          fprintf(out,"\nDRIVER: Reassigning domain decomposition\n");
        }
        p.cmd("setAtomsNlocal",&dd_nlocal);
        p.cmd("setAtomsGatindex",&dd_gatindex[0]);
      }
    }

    int plumedStopCondition=0;
    if(!noatoms) {
      if(use_molfile) {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
        if(pbc_cli_given==false) {
          if(ts_in.A>0.0) { // this is negative if molfile does not provide box
            // info on the cell: convert using pbcset.tcl from pbctools in vmd distribution
            real cosBC=cos(real(ts_in.alpha)*pi/180.);
            //double sinBC=sin(ts_in.alpha*pi/180.);
            real cosAC=cos(real(ts_in.beta)*pi/180.);
            real cosAB=cos(real(ts_in.gamma)*pi/180.);
            real sinAB=sin(real(ts_in.gamma)*pi/180.);
            real Ax=real(ts_in.A);
            real Bx=real(ts_in.B)*cosAB;
            real By=real(ts_in.B)*sinAB;
            real Cx=real(ts_in.C)*cosAC;
            real Cy=(real(ts_in.C)*real(ts_in.B)*cosBC-Cx*Bx)/By;
            real Cz=sqrt(real(ts_in.C)*real(ts_in.C)-Cx*Cx-Cy*Cy);
            cell[0]=Ax/10.; cell[1]=0.; cell[2]=0.;
            cell[3]=Bx/10.; cell[4]=By/10.; cell[5]=0.;
            cell[6]=Cx/10.; cell[7]=Cy/10.; cell[8]=Cz/10.;
          } else {
            cell[0]=0.0; cell[1]=0.0; cell[2]=0.0;
            cell[3]=0.0; cell[4]=0.0; cell[5]=0.0;
            cell[6]=0.0; cell[7]=0.0; cell[8]=0.0;
          }
        } else {
          for(unsigned i=0; i<9; i++)cell[i]=pbc_cli_box[i];
        }
        // info on coords
        // the order is xyzxyz...
        for(int i=0; i<3*natoms; i++) {
          coordinates[i]=real(ts_in.coords[i])/real(10.); //convert to nm
          //cerr<<"COOR "<<coordinates[i]<<endl;
        }
#endif
      } else if(trajectory_fmt=="xdr-xtc" || trajectory_fmt=="xdr-trr") {
#ifdef __PLUMED_HAS_XDRFILE
        int localstep;
        float time;
        matrix box;
        std::unique_ptr<rvec[]> pos(new rvec[natoms]);
        float prec,lambda;
        int ret=exdrOK;
        if(trajectory_fmt=="xdr-xtc") ret=read_xtc(xd,natoms,&localstep,&time,box,pos.get(),&prec);
        if(trajectory_fmt=="xdr-trr") ret=read_trr(xd,natoms,&localstep,&time,&lambda,box,pos.get(),NULL,NULL);
        if(stride==0) step=localstep;
        if(ret==exdrENDOFFILE) break;
        if(ret!=exdrOK) break;
        for(unsigned i=0; i<3; i++) for(unsigned j=0; j<3; j++) cell[3*i+j]=box[i][j];
        for(int i=0; i<natoms; i++) for(unsigned j=0; j<3; j++)
            coordinates[3*i+j]=real(pos[i][j]);
#endif
      } else {
        if(trajectory_fmt=="xyz") {
          if(!Tools::getline(fp,line)) error("premature end of trajectory file");

          std::vector<double> celld(9,0.0);
          if(pbc_cli_given==false) {
            std::vector<std::string> words;
            words=Tools::getWords(line);
            if(words.size()==3) {
              sscanf(line.c_str(),"%100lf %100lf %100lf",&celld[0],&celld[4],&celld[8]);
            } else if(words.size()==9) {
              sscanf(line.c_str(),"%100lf %100lf %100lf %100lf %100lf %100lf %100lf %100lf %100lf",
                     &celld[0], &celld[1], &celld[2],
                     &celld[3], &celld[4], &celld[5],
                     &celld[6], &celld[7], &celld[8]);
            } else error("needed box in second line of xyz file");
          } else {			// from command line
            celld=pbc_cli_box;
          }
          for(unsigned i=0; i<9; i++)cell[i]=real(celld[i]);
        }
        if(trajectory_fmt=="dlp4") {
          std::vector<double> celld(9,0.0);
          if(pbc_cli_given==false) {
            if(!Tools::getline(fp,line)) error("error reading vector a of cell");
            sscanf(line.c_str(),"%lf %lf %lf",&celld[0],&celld[1],&celld[2]);
            if(!Tools::getline(fp,line)) error("error reading vector b of cell");
            sscanf(line.c_str(),"%lf %lf %lf",&celld[3],&celld[4],&celld[5]);
            if(!Tools::getline(fp,line)) error("error reading vector c of cell");
            sscanf(line.c_str(),"%lf %lf %lf",&celld[6],&celld[7],&celld[8]);
          } else {
            celld=pbc_cli_box;
          }
          for(auto i=0; i<9; i++)cell[i]=real(celld[i])*0.1;
        }
        int ddist=0;
        // Read coordinates
        for(int i=0; i<natoms; i++) {
          bool ok=Tools::getline(fp,line);
          if(!ok) error("premature end of trajectory file");
          double cc[3];
          if(trajectory_fmt=="xyz") {
            char dummy[1000];
            int ret=std::sscanf(line.c_str(),"%999s %100lf %100lf %100lf",dummy,&cc[0],&cc[1],&cc[2]);
            if(ret!=4) error("cannot read line"+line);
          } else if(trajectory_fmt=="gro") {
            // do the gromacs way
            if(!i) {
              //
              // calculate the distance between dots (as in gromacs gmxlib/confio.c, routine get_w_conf )
              //
              const char      *p1, *p2, *p3;
              p1 = strchr(line.c_str(), '.');
              if (p1 == NULL) error("seems there are no coordinates in the gro file");
              p2 = strchr(&p1[1], '.');
              if (p2 == NULL) error("seems there is only one coordinates in the gro file");
              ddist = p2 - p1;
              p3 = strchr(&p2[1], '.');
              if (p3 == NULL)error("seems there are only two coordinates in the gro file");
              if (p3 - p2 != ddist)error("not uniform spacing in fields in the gro file");
            }
            Tools::convert(line.substr(20,ddist),cc[0]);
            Tools::convert(line.substr(20+ddist,ddist),cc[1]);
            Tools::convert(line.substr(20+ddist+ddist,ddist),cc[2]);
          } else if(trajectory_fmt=="dlp4") {
            char dummy[9];
            int idummy;
            double m,c;
            sscanf(line.c_str(),"%8s %d %lf %lf",dummy,&idummy,&m,&c);
            masses[i]=real(m);
            charges[i]=real(c);
            if(!Tools::getline(fp,line)) error("error reading coordinates");
            sscanf(line.c_str(),"%lf %lf %lf",&cc[0],&cc[1],&cc[2]);
            cc[0]*=0.1;
            cc[1]*=0.1;
            cc[2]*=0.1;
            if(lvl>0) {
              if(!Tools::getline(fp,line)) error("error skipping velocities");
            }
            if(lvl>1) {
              if(!Tools::getline(fp,line)) error("error skipping forces");
            }
          } else plumed_error();
          if(!debug_pd || ( i>=pd_start && i<pd_start+pd_nlocal) ) {
            coordinates[3*i]=real(cc[0]);
            coordinates[3*i+1]=real(cc[1]);
            coordinates[3*i+2]=real(cc[2]);
          }
        }
        if(trajectory_fmt=="gro") {
          if(!Tools::getline(fp,line)) error("premature end of trajectory file");
          std::vector<string> words=Tools::getWords(line);
          if(words.size()<3) error("cannot understand box format");
          Tools::convert(words[0],cell[0]);
          Tools::convert(words[1],cell[4]);
          Tools::convert(words[2],cell[8]);
          if(words.size()>3) Tools::convert(words[3],cell[1]);
          if(words.size()>4) Tools::convert(words[4],cell[2]);
          if(words.size()>5) Tools::convert(words[5],cell[3]);
          if(words.size()>6) Tools::convert(words[6],cell[5]);
          if(words.size()>7) Tools::convert(words[7],cell[6]);
          if(words.size()>8) Tools::convert(words[8],cell[7]);
        }

      }

      p.cmd("setStepLong",&step);
      p.cmd("setStopFlag",&plumedStopCondition);

      if(debug_dd) {
        for(int i=0; i<dd_nlocal; ++i) {
          int kk=dd_gatindex[i];
          dd_coordinates[3*i+0]=coordinates[3*kk+0];
          dd_coordinates[3*i+1]=coordinates[3*kk+1];
          dd_coordinates[3*i+2]=coordinates[3*kk+2];
        }
        p.cmd("setForces",&dd_forces[0]);
        p.cmd("setPositions",&dd_coordinates[0]);
        p.cmd("setMasses",&dd_masses[0]);
        p.cmd("setCharges",&dd_charges[0]);
      } else {
// this is required to avoid troubles when the last domain
// contains zero atoms
// Basically, for empty domains we pass null pointers
#define fix_pd(xx) (pd_nlocal!=0?&xx:NULL)
        p.cmd("setForces",fix_pd(forces[3*pd_start]));
        p.cmd("setPositions",fix_pd(coordinates[3*pd_start]));
        p.cmd("setMasses",fix_pd(masses[pd_start]));
        p.cmd("setCharges",fix_pd(charges[pd_start]));
      }
      p.cmd("setBox",&cell[0]);
      p.cmd("setVirial",&virial[0]);
    } else {
      p.cmd("setStepLong",&step);
      p.cmd("setStopFlag",&plumedStopCondition);
    }
    p.cmd("calc");
    if(debugforces.length()>0) {
      virial.assign(9,real(0.0));
      forces.assign(3*natoms,real(0.0));
      p.cmd("prepareCalc");
      p.cmd("performCalcNoUpdate");
    }

// this is necessary as only processor zero is adding to the virial:
    intracomm.Bcast(virial,0);
    if(debug_pd) intracomm.Sum(forces);
    if(debug_dd) {
      for(int i=0; i<dd_nlocal; i++) {
        forces[3*dd_gatindex[i]+0]=dd_forces[3*i+0];
        forces[3*dd_gatindex[i]+1]=dd_forces[3*i+1];
        forces[3*dd_gatindex[i]+2]=dd_forces[3*i+2];
      }
      dd_forces.assign(3*natoms,0.0);
      intracomm.Sum(forces);
    }
    if(debug_grex &&step%grex_stride==0) {
      p.cmd("GREX savePositions");
      if(intracomm.Get_rank()>0) {
        p.cmd("GREX prepare");
      } else {
        int r=intercomm.Get_rank();
        int n=intercomm.Get_size();
        int partner=r+(2*((r+step/grex_stride)%2))-1;
        if(partner<0)partner=0;
        if(partner>=n) partner=n-1;
        p.cmd("GREX setPartner",&partner);
        p.cmd("GREX calculate");
        p.cmd("GREX shareAllDeltaBias");
        for(int i=0; i<n; i++) {
          string s; Tools::convert(i,s);
          real a=std::numeric_limits<real>::quiet_NaN(); s="GREX getDeltaBias "+s; p.cmd(s.c_str(),&a);
          if(grex_log) fprintf(grex_log," %f",a);
        }
        if(grex_log) fprintf(grex_log,"\n");
      }
    }


    if(fp_forces) {
      fprintf(fp_forces,"%d\n",natoms);
      string fmtv=dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+"\n";
      string fmt=dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+"\n";
      if(dumpfullvirial) {
        fprintf(fp_forces,fmtv.c_str(),virial[0],virial[1],virial[2],virial[3],virial[4],virial[5],virial[6],virial[7],virial[8]);
      } else {
        fprintf(fp_forces,fmt.c_str(),virial[0],virial[4],virial[8]);
      }
      fmt="X "+fmt;
      for(int i=0; i<natoms; i++)
        fprintf(fp_forces,fmt.c_str(),forces[3*i],forces[3*i+1],forces[3*i+2]);
    }
    if(debugforces.length()>0) {
      // Now call the routine to work out the derivatives numerically
      numder.assign(3*natoms+9,real(0.0)); real base=0;
      p.cmd("getBias",&base);
      if( fabs(base)<epsilon ) printf("WARNING: bias for configuration appears to be zero so debugging forces is trivial");
      evaluateNumericalDerivatives( step, p, coordinates, masses, charges, cell, base, numder );

      // And output everything to a file
      fp_dforces.fmtField(" " + dumpforcesFmt);
      for(int i=0; i<3*natoms; ++i) {
        fp_dforces.printField("parameter",(int)i);
        fp_dforces.printField("analytical",forces[i]);
        fp_dforces.printField("numerical",-numder[i]);
        fp_dforces.printField();
      }
      // And print the virial
      for(int i=0; i<9; ++i) {
        fp_dforces.printField("parameter",3*natoms+i );
        fp_dforces.printField("analytical",virial[i] );
        fp_dforces.printField("numerical",-numder[3*natoms+i]);
        fp_dforces.printField();
      }
    }

    if(plumedStopCondition) break;

    step+=stride;
  }
  if(!parseOnly) p.cmd("runFinalJobs");

  if(fp_forces) fclose(fp_forces);
  if(debugforces.length()>0) fp_dforces.close();
  if(fp && fp!=in)fclose(fp);
#ifdef __PLUMED_HAS_XDRFILE
  if(xd) xdrfile_close(xd);
#endif
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
  if(h_in) api->close_file_read(h_in);
#endif
  if(grex_log) fclose(grex_log);

  return 0;
}

template<typename real>
void Driver<real>::evaluateNumericalDerivatives( const long int& step, PlumedMain& p, const std::vector<real>& coordinates,
    const std::vector<real>& masses, const std::vector<real>& charges,
    std::vector<real>& cell, const double& base, std::vector<real>& numder ) {

  int natoms = coordinates.size() / 3; real delta = sqrt(epsilon);
  std::vector<Vector> pos(natoms); real bias=0;
  std::vector<real> fake_forces( 3*natoms ), fake_virial(9);
  for(int i=0; i<natoms; ++i) {
    for(unsigned j=0; j<3; ++j) pos[i][j]=coordinates[3*i+j];
  }

  for(int i=0; i<natoms; ++i) {
    for(unsigned j=0; j<3; ++j) {
      pos[i][j]=pos[i][j]+delta;
      p.cmd("setStepLong",&step);
      p.cmd("setPositions",&pos[0][0]);
      p.cmd("setForces",&fake_forces[0]);
      p.cmd("setMasses",&masses[0]);
      p.cmd("setCharges",&charges[0]);
      p.cmd("setBox",&cell[0]);
      p.cmd("setVirial",&fake_virial[0]);
      p.cmd("prepareCalc");
      p.cmd("performCalcNoUpdate");
      p.cmd("getBias",&bias);
      pos[i][j]=coordinates[3*i+j];
      numder[3*i+j] = (bias - base) / delta;
    }
  }

  // Create the cell
  Tensor box( cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8] );
  // And the virial
  Pbc pbc; pbc.setBox( box ); Tensor nvirial;
  for(unsigned i=0; i<3; i++) for(unsigned k=0; k<3; k++) {
      double arg0=box(i,k);
      for(int j=0; j<natoms; ++j) pos[j]=pbc.realToScaled( pos[j] );
      cell[3*i+k]=box(i,k)=box(i,k)+delta; pbc.setBox(box);
      for(int j=0; j<natoms; j++) pos[j]=pbc.scaledToReal( pos[j] );
      p.cmd("setStepLong",&step);
      p.cmd("setPositions",&pos[0][0]);
      p.cmd("setForces",&fake_forces[0]);
      p.cmd("setMasses",&masses[0]);
      p.cmd("setCharges",&charges[0]);
      p.cmd("setBox",&cell[0]);
      p.cmd("setVirial",&fake_virial[0]);
      p.cmd("prepareCalc");
      p.cmd("performCalcNoUpdate");
      p.cmd("getBias",&bias);
      cell[3*i+k]=box(i,k)=arg0; pbc.setBox(box);
      for(int j=0; j<natoms; j++) for(unsigned n=0; n<3; ++n) pos[j][n]=coordinates[3*j+n];
      nvirial(i,k) = ( bias - base ) / delta;
    }
  nvirial=-matmul(box.transpose(),nvirial);
  for(unsigned i=0; i<3; i++) for(unsigned k=0; k<3; k++)  numder[3*natoms+3*i+k] = nvirial(i,k);

}

}
}
