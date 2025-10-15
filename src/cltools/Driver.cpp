/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/CLToolRegister.h"
#include "tools/Tools.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithVirtualAtom.h"
#include "core/ActionShortcut.h"
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
#include "tools/TrajectoryParser.h"

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS driver-float
/*
A tool that does the same as driver, but using single precision reals.

We recommend always using [driver](driver.md) to postprocess your trajectories. However,
if you want to test what PLUMED does when linked from a single precision code you can use
this version of driver. The documentation for this command line tool is identical to that for
[driver](driver.md).

## Examples

```plumed
plumed driver-float --plumed plumed.dat --ixyz trajectory.xyz
```

For more examples see the documentation for [driver](driver.md).

*/
//+ENDPLUMEDOC
//


//+PLUMEDOC TOOLS driver
/*
driver is a tool that allows one to to use plumed to post-process an existing trajectory.

The input to driver is specified using the command line arguments described below.

In addition, you can use the special [READ](READ.md) command inside your plumed input
to read in colvar files that were generated during your MD simulation.  The values
read in from the file can then be treated like calculated colvars.

!!! warning "masses and chages"

    Notice that by default the driver has no knowledge about the masses and charges
    of your atoms! Thus, if you want to compute quantities that depend on charges (e.g. [DHENERGY](DHENERGY.md))
    or masses (e.g. [COM](COM.md)) you need to pass masses and charges to driver.
    You can do this by either using --pdb option or the --mc option. The latter
    will read a file produced by [DUMPMASSCHARGE](DUMPMASSCHARGE.md).


## Examples

The following command tells plumed to post process the trajectory contained in `trajectory.xyz`
 by performing the actions described in the input file `plumed.dat`.  If an action that takes the
stride keyword is given a stride equal to $n$ then it will be performed only on every $n$th
frame in the trajectory file.

```plumed
plumed driver --plumed plumed.dat --ixyz trajectory.xyz
```

Notice that `xyz` files are expected to be in the internal PLUMED unit of nm.
You can change this behavior by using the `--length-units` option as shown below

```plumed
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --length-units A
```

The strings accepted by the `--length-units` options are the same ones accepted by the [UNITS](UNITS.md) action.
Other file formats typically have their default coordinates (e.g., `gro` files are always in nm)
and it thus should not be necessary to use the `--length-units` option. Additionally,
consider that the units used by the `driver` might be different by the units used in the PLUMED input
file `plumed.dat`. For instance consider the command:

```plumed
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --length-units A
```

which uses the following `plumed.dat` file

```plumed
# no explicit UNITS action here
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=colvar
```

In this case, the driver reads the `xyz` file assuming it to contain coordinates in Angstroms.
However, the resulting `colvar` file contains a distance expressed in nm.

The following command tells plumed to post process the trajectory contained in trajectory.xyz.
 by performing the actions described in the input file plumed.dat.

```plumed
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --trajectory-stride 100 --timestep 0.001
```

Here though
`--trajectory-stride` is set equal to the frequency with which frames were output during the trajectory
and the `--timestep` is equal to the simulation timestep.  As such the `STRIDE` parameters in the `plumed.dat`
files are referred to the original timestep and any files output resemble those that would have been generated
had we run the calculation we are running with driver when the MD simulation was running.

PLUMED can read xyz files (in PLUMED units) and gro files (in nm). In addition,
PLUMED includes support for a
subset of the trajectory file formats supported by VMD, e.g. xtc and dcd:

```plumed
plumed driver --plumed plumed.dat --pdb diala.pdb --mf_xtc traj.xtc --trajectory-stride 100 --timestep 0.001
```

where `--mf_` prefixes the extension of one of the accepted molfile plugin format.
If PLUMED has been installed with full molfile support, other formats will be available.
Just type `plumed driver --help` to see which file formats you can use.

Molfile plugin requires the periodic cell to be triangular (i.e. the first vector oriented along x and
second vector in xy plane). This is true for many MD codes. However, it could be false
if you rotate the coordinates in your trajectory before reading them in the driver.
Also notice that some formats (e.g. amber crd) do not specify the total number of atoms. In this case you can use
the `--natoms` option as shown below

```plumed
plumed driver --plumed plumed.dat --mf_crd trajectory.crd --natoms 128
```

Check the available molfile plugins and limitations at [this link](http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/).

You can also use the xdrfile implementation of xtc and trr with plumed driver. If you want to do so you just
download and install properly the xdrfile library (see [this link](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library)).
If the xdrfile library is installed properly the PLUMED configure script should be able to
detect it and enable it.
Notice that the xdrfile implementation of xtc and trr
is more robust than the molfile one, since it provides support for generic cell shapes.
In addition, if you install xdrfile you can then use the [DUMPATOMS](DUMPATOMS.md) command to write compressed xtc files.

## Multiple replicas

When PLUMED is compiled with MPI support, you can emulate a multi-simulation setup with `driver` by providing the `--multi`
option with the appropriate number of ranks. This allows you to use the ref special-replica-syntax that is discussed [here](parsing.md) to analyze multiple
trajectories (see [this tutorial](https://www.plumed-tutorials.org/lessons/21/005/data/NAVIGATION.html)). PLUMED will also automatically append a numbered suffix to output files
(e.g. `COLVAR.0`, `COLVAR.1`, …) as discussed [here](Files.md). Similarly, each replica will search for the corresponding suffixed input file (e.g. `traj.0.xtc`, …)
or default to the unsuffixed one.

*/
//+ENDPLUMEDOC
//

template<typename real>
class Driver : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit Driver(const CLToolOptions& co );
  int main(FILE* in,FILE*out,Communicator& pc) override;
  void evaluateNumericalDerivatives( const long long int& step, PlumedMain& p, const std::vector<real>& coordinates,
                                     const std::vector<real>& masses, const std::vector<real>& charges,
                                     std::vector<real>& cell, const double& base, std::vector<real>& numder );
  std::string description()const override;
};

template<typename real>
void Driver<real>::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.isDriver();
  keys.addFlag("--help-debug",false,"print special options that can be used to create regtests");
  keys.add("compulsory","--plumed","plumed.dat","specify the name of the plumed input file");
  keys.add("compulsory","--timestep","1.0","the timestep that was used in the calculation that produced this trajectory in picoseconds");
  keys.add("compulsory","--trajectory-stride","1","the frequency with which frames were output to this trajectory during the simulation"
           " (0 means that the number of the step is read from the trajectory file,"
           " currently working only for xtc/trr files read with --ixtc/--trr)"
          );
  keys.add("compulsory","--multi","0","set number of replicas for multi environment (needs MPI)");
  keys.addFlag("--noatoms",false,"don't read in a trajectory.  Just use colvar files as specified in plumed.dat");
  keys.addFlag("--parse-only",false,"read the plumed input file and stop");
  keys.addFlag("--restart",false,"makes driver behave as if restarting");
  TrajectoryParser::registerKeywords(keys);
  keys.add("optional","--shortcut-ofile","the name of the file to output info on the way shortcuts have been expanded.  If there are no shortcuts in your input file nothing is output");
  keys.add("optional","--valuedict-ofile","output a dictionary giving information about each value in the input file");
  keys.add("optional","--length-units","units for length, either as a string or a number");
  keys.add("optional","--mass-units","units for mass in pdb and mc file, either as a string or a number");
  keys.add("optional","--charge-units","units for charge in pdb and mc file, either as a string or a number");
  keys.add("optional","--kt","set kT it will not be necessary to specify temperature in input file");
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
}
template<typename real>
Driver<real>::Driver(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}
template<typename real>
std::string Driver<real>::description()const {
  return "analyze trajectories with plumed";
}

template<typename real>
int Driver<real>::main(FILE* in,FILE*out,Communicator& pc) {

  Units units;
  PDB pdb;

// Parse everything
  bool printhelpdebug;
  parseFlag("--help-debug",printhelpdebug);
  if( printhelpdebug ) {
    std::fprintf(out,"%s",
                 "Additional options for debug (only to be used in regtest):\n"
                 "  [--debug-float yes]     : turns on the single precision version (to check float interface)\n"
                 "  [--debug-dd yes]        : use a fake domain decomposition\n"
                 "  [--debug-pd yes]        : use a fake particle decomposition\n"
                );
    return 0;
  }
  // Are we reading trajectory data
  bool noatoms;
  parseFlag("--noatoms",noatoms);
  bool parseOnly;
  parseFlag("--parse-only",parseOnly);
  std::string full_outputfile;
  parse("--shortcut-ofile",full_outputfile);
  std::string valuedict_file;
  parse("--valuedict-ofile",valuedict_file);
  bool restart;
  parseFlag("--restart",restart);

  std::string fakein;
  bool debug_float=false;
  fakein="";
  if(parse("--debug-float",fakein)) {
    if(fakein=="yes") {
      debug_float=true;
    } else if(fakein=="no") {
      debug_float=false;
    } else {
      error("--debug-float should have argument yes or no");
    }
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
    if(fakein=="yes") {
      debug_pd=true;
    } else if(fakein=="no") {
      debug_pd=false;
    } else {
      error("--debug-pd should have argument yes or no");
    }
  }
  if(debug_pd) {
    std::fprintf(out,"DEBUGGING PARTICLE DECOMPOSITION\n");
  }

  bool debug_dd=false;
  fakein="";
  if(parse("--debug-dd",fakein)) {
    if(fakein=="yes") {
      debug_dd=true;
    } else if(fakein=="no") {
      debug_dd=false;
    } else {
      error("--debug-dd should have argument yes or no");
    }
  }
  if(debug_dd) {
    std::fprintf(out,"DEBUGGING DOMAIN DECOMPOSITION\n");
  }

  if( debug_pd || debug_dd ) {
    if(noatoms) {
      error("cannot debug without atoms");
    }
  }

// set up for multi replica driver:
  int multi=0;
  parse("--multi",multi);
  Communicator intracomm;
  Communicator intercomm;
  if(multi) {
    int ntot=pc.Get_size();
    int nintra=ntot/multi;
    if(multi*nintra!=ntot) {
      error("invalid number of processes for multi environment");
    }
    pc.Split(pc.Get_rank()/nintra,pc.Get_rank(),intracomm);
    pc.Split(pc.Get_rank()%nintra,pc.Get_rank(),intercomm);
  } else {
    intracomm.Set_comm(pc.Get_comm());
  }

// set up for debug replica exchange:
  bool debug_grex=parse("--debug-grex",fakein);
  int  grex_stride=0;
  FILE*grex_log=NULL;
// call fclose when fp goes out of scope
  auto deleter=[](auto f) {
    if(f) {
      std::fclose(f);
    }
  };
  std::unique_ptr<FILE,decltype(deleter)> grex_log_deleter(grex_log,deleter);

  if(debug_grex) {
    if(noatoms) {
      error("must have atoms to debug_grex");
    }
    if(multi<2) {
      error("--debug_grex needs --multi with at least two replicas");
    }
    Tools::convert(fakein,grex_stride);
    std::string n;
    Tools::convert(intercomm.Get_rank(),n);
    std::string file;
    parse("--debug-grex-log",file);
    if(file.length()>0) {
      file+="."+n;
      grex_log=std::fopen(file.c_str(),"w");
      grex_log_deleter.reset(grex_log);
    }
  }

// Read the plumed input file name
  std::string plumedFile;
  parse("--plumed",plumedFile);
// the timestep
  double t;
  parse("--timestep",t);
  real timestep=real(t);
// the stride
  unsigned stride;
  parse("--trajectory-stride",stride);
// are we writing forces
  std::string dumpforces(""), debugforces(""), dumpforcesFmt("%f");;
  bool dumpfullvirial=false;
  if(!noatoms) {
    parse("--dump-forces",dumpforces);
    parse("--debug-forces",debugforces);
  }
  if(dumpforces!="" || debugforces!="" ) {
    parse("--dump-forces-fmt",dumpforcesFmt);
  }
  if(dumpforces!="") {
    parseFlag("--dump-full-virial",dumpfullvirial);
  }
  if( debugforces!="" && (debug_dd || debug_pd) ) {
    error("cannot debug forces and domain/particle decomposition at same time");
  }
  if( debugforces!="" && sizeof(real)!=sizeof(double) ) {
    error("cannot debug forces in single precision mode");
  }

  real kt=-1.0;
  parse("--kt",kt);
  std::string trajectory_fmt;

  bool use_molfile=false;

  // Read in an xyz file
  std::string trajectoryFile(""), pdbfile(""), mcfile("");
  bool pbc_cli_given=false;
  std::vector<double> pbc_cli_box(9,0.0);
  int command_line_natoms=-1;

  TrajectoryParser parser;
  if(!noatoms) {
    int nn=0;
    trajectory_fmt="";
    for (const auto & trj_type : TrajectoryParser::trajectoryOptions()) {
      std::string tmp;
      parse("--i"+trj_type, tmp);
      if (tmp.length()>0) {
        trajectory_fmt=trj_type;
        ++nn;
        trajectoryFile=tmp;
      }
    }
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
    {
      auto plugins_names=TrajectoryParser::getMolfilePluginsnames() ;
      for(unsigned i=0; i<plugins_names.size(); i++) {
        std::string molfile_key="--mf_"+plugins_names[i];
        std::string traj_molfile;
        parse(molfile_key,traj_molfile);
        if(traj_molfile.length()>0) {
          ++nn;
          std::fprintf(out,"\nDRIVER: Found molfile format trajectory %s with name %s\n",plugins_names[i].c_str(),traj_molfile.c_str());
          trajectoryFile=traj_molfile;
          trajectory_fmt=plugins_names[i];
          use_molfile=true;
        }
      }
    }
#endif
    {

      if(nn>1) {
        std::fprintf(stderr,"ERROR: cannot provide more than one trajectory file\n");
        return 1;
      }
    }

    if(trajectoryFile.length()==0&&!parseOnly) {
      std::fprintf(stderr,"ERROR: missing trajectory data\n");
      return 1;
    }
    std::string lengthUnits("");
    parse("--length-units",lengthUnits);
    if(lengthUnits.length()>0) {
      units.setLength(lengthUnits);
    }
    std::string chargeUnits("");
    parse("--charge-units",chargeUnits);
    if(chargeUnits.length()>0) {
      units.setCharge(chargeUnits);
    }
    std::string massUnits("");
    parse("--mass-units",massUnits);
    if(massUnits.length()>0) {
      units.setMass(massUnits);
    }

    parse("--pdb",pdbfile);
    if(pdbfile.length()>0) {
      bool check=pdb.read(pdbfile,false,1.0);
      if(!check) {
        error("error reading pdb file");
      }
    }

    parse("--mc",mcfile);

    std::string pbc_cli_list;
    parse("--box",pbc_cli_list);
    if(pbc_cli_list.length()>0) {
      pbc_cli_given=true;
      std::vector<std::string> words=Tools::getWords(pbc_cli_list,",");
      if(words.size()==3) {
        for(int i=0; i<3; i++) {
          std::sscanf(words[i].c_str(),"%100lf",&(pbc_cli_box[4*i]));
        }
      } else if(words.size()==9) {
        for(int i=0; i<9; i++) {
          std::sscanf(words[i].c_str(),"%100lf",&(pbc_cli_box[i]));
        }
      } else {
        std::string msg="ERROR: cannot parse command-line box "+pbc_cli_list;
        std::fprintf(stderr,"%s\n",msg.c_str());
        return 1;
      }

    }

    parse("--natoms",command_line_natoms);

  } // if(!noatoms)

  if(debug_dd && debug_pd) {
    error("cannot use debug-dd and debug-pd at the same time");
  }
  if(debug_pd || debug_dd) {
    if( !Communicator::initialized() ) {
      error("needs mpi for debug-pd");
    }
  }

  PlumedMain p;
  if( parseOnly ) {
    p.activateParseOnlyMode();
  }
  p.cmd("setRealPrecision",(int)sizeof(real));
  int checknatoms=-1;
  long long int step=0;
  parse("--initial-step",step);

  if(restart) {
    p.cmd("setRestart",1);
  }

  if(Communicator::initialized()) {
    if(multi) {
      if(intracomm.Get_rank()==0) {
        p.cmd("GREX setMPIIntercomm",&intercomm.Get_comm());
      }
      p.cmd("GREX setMPIIntracomm",&intracomm.Get_comm());
      p.cmd("GREX init");
    }
    p.cmd("setMPIComm",&intracomm.Get_comm());
  }
  p.cmd("setLog",out);
  p.cmd("setMDLengthUnits",units.getLength());
  p.cmd("setMDChargeUnits",units.getCharge());
  p.cmd("setMDMassUnits",units.getMass());
  p.cmd("setMDEngine","driver");
  p.cmd("setTimestep",timestep);
  if( !parseOnly || full_outputfile.length()==0 ) {
    p.cmd("setPlumedDat",plumedFile.c_str());
  }

  int natoms=0;

  if(parseOnly) {
    if(command_line_natoms<0) {
      error("--parseOnly requires setting the number of atoms with --natoms");
    }
    natoms=command_line_natoms;
  }


  FILE* fp_forces=NULL;
  OFile fp_dforces;

  std::unique_ptr<FILE,decltype(deleter)> fp_forces_deleter(fp_forces,deleter);

  if(!noatoms&&!parseOnly) {
    if (trajectoryFile=="-") {
      parser.init(trajectory_fmt, in);
    } else {
      if(multi) {
        std::string n;
        Tools::convert(intercomm.Get_rank(),n);
        std::string testfile=FileBase::appendSuffix(trajectoryFile,"."+n);
        FILE* tmp_fp=std::fopen(testfile.c_str(),"r");
        // no exceptions here
        if(tmp_fp) {
          std::fclose(tmp_fp);
          trajectoryFile=testfile;
        }
      }
      parser.init(trajectory_fmt,trajectoryFile,use_molfile,command_line_natoms);
      natoms = parser.nOfAtoms();
    }
    if(dumpforces.length()>0) {
      if(Communicator::initialized() && pc.Get_size()>1) {
        std::string n;
        Tools::convert(pc.Get_rank(),n);
        dumpforces+="."+n;
      }
      fp_forces=std::fopen(dumpforces.c_str(),"w");
      fp_forces_deleter.reset(fp_forces);
    }
    if(debugforces.length()>0) {
      if(Communicator::initialized() && pc.Get_size()>1) {
        std::string n;
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
  int pd_nlocal=0;
  int pd_start=0;
// variables to test random decomposition (=domain decomposition)
  std::vector<int>  dd_gatindex;
  std::vector<int>  dd_g2l;
  std::vector<real> dd_masses;
  std::vector<real> dd_charges;
  std::vector<real> dd_forces;
  std::vector<real> dd_coordinates;
  int dd_nlocal=0;
// random stream to choose decompositions
  Random rnd;


  bool lstep=true;
  while(true) {
    bool first_step=false;
    if(!noatoms&&!parseOnly) {

      real timeStep=-1.0;

      auto errormessage=parser.readHeader(step,timeStep);

      if (errormessage) {
        if (*errormessage =="EOF") {
          break;
        } else {
          error(*errormessage);
        }
      }
      if (lstep && timeStep>0.0) {
        p.cmd("setTimestep",real(timeStep));
        lstep = false;
      }
      natoms=parser.nOfAtoms();
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
          if( index>=unsigned(natoms) ) {
            error("atom index in pdb exceeds the number of atoms in trajectory");
          }
          masses[index]=pdb.getOccupancy()[i];
          charges[index]=pdb.getBeta()[i];
        }
      }
      if(mcfile.length()>0) {
        IFile ifile;
        ifile.open(mcfile);
        int index;
        double mass;
        double charge;
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
        p.cmd("setKbT",kt);
      }
      checknatoms=natoms;
      p.cmd("setNatoms",natoms);
      p.cmd("init");
      // Check if we have been asked to output the long version of the input and if there are shortcuts
      if( parseOnly && full_outputfile.length()>0 ) {

        // Read in the plumed input file and store what is in there
        std::map<std::string,std::vector<std::string> > data;
        IFile ifile;
        ifile.open(plumedFile);
        std::vector<std::string> words;
        while( Tools::getParsedLine(ifile,words) && !p.getEndPlumed() ) {
          p.readInputWords(words,false);
          Action* aa=p.getActionSet()[p.getActionSet().size()-1].get();
          plumed_assert(aa); // needed for following calls, see #1046
          ActionWithValue* av=aa->castToActionWithValue();
          if( av && aa->getDefaultString().length()>0 ) {
            std::vector<std::string> def;
            def.push_back( "defaults " + aa->getDefaultString() );
            data[ aa->getLabel() ] = def;
          }
          ActionShortcut* as=aa->castToActionShortcut();
          if( as ) {
            if( aa->getDefaultString().length()>0 ) {
              std::vector<std::string> def;
              def.push_back( "defaults " + aa->getDefaultString() );
              data[ as->getShortcutLabel() ] = def;
            }
            const auto & shortcut_commands = as->getSavedInputLines();
            if( shortcut_commands.size()==0 ) {
              continue;
            }
            if( data.find( as->getShortcutLabel() )!=data.end() ) {
              for(unsigned i=0; i<shortcut_commands.size(); ++i) {
                data[ as->getShortcutLabel() ].push_back( shortcut_commands[i] );
              }
            } else {
              data[ as->getShortcutLabel() ] = shortcut_commands;
            }
          }
        }
        ifile.close();
        // Only output the full version of the input file if there are shortcuts
        if( data.size()>0 && intracomm.Get_rank()==0 && intercomm.Get_rank()==0 ) {
          OFile long_file;
          long_file.open( full_outputfile );
          long_file.printf("{\n");
          bool firstpass=true;
          for(auto& x : data ) {
            if( !firstpass ) {
              long_file.printf("   },\n");
            }
            long_file.printf("   \"%s\" : {\n", x.first.c_str() );
            plumed_assert( x.second.size()>0 );
            unsigned sstart=0;
            if( x.second[0].find("defaults")!=std::string::npos ) {
              sstart=1;
              long_file.printf("      \"defaults\" : \"%s\"", x.second[0].substr( 9 ).c_str() );
              if( x.second.size()>1 ) {
                long_file.printf(",\n");
              } else {
                long_file.printf("\n");
              }
            }
            if( x.second.size()>sstart ) {
              {
                std::string cmnd=Tools::convertRegexForJson(x.second[sstart]);
                long_file.printf("      \"expansion\" : \"%s", cmnd.c_str() );
              }
              for(unsigned j=sstart+1; j<x.second.size(); ++j) {
                //reinitializing it to RVO the string
                std::string cmnd=Tools::convertRegexForJson(x.second[j]);
                long_file.printf("\\n%s", cmnd.c_str() );
              }
              long_file.printf("\"\n");
            }
            firstpass=false;
          }
          long_file.printf("   }\n}\n");
          long_file.close();
        }
      }
      if( valuedict_file.length()>0 && intracomm.Get_rank()==0 && intercomm.Get_rank()==0 ) {
        OFile valuefile;
        valuefile.open( valuedict_file );
        valuefile.printf("{\n");
        bool firsta=true;
        for(const auto & pp : p.getActionSet()) {
          if( pp.get()->getName()=="CENTER" ) {
            continue ;
          }
          ActionWithVirtualAtom* avv=dynamic_cast<ActionWithVirtualAtom*>(pp.get());
          if( avv ||  pp.get()->getName()=="GROUP" || pp.get()->getName()=="DENSITY" ) {
            Action* actionBase=pp.get();
            if( firsta ) {
              valuefile.printf("  \"%s\" : {\n    \"action\" : \"%s\"",
                               actionBase->getLabel().c_str(),
                               actionBase->getName().c_str() );
              firsta=false;
            } else {
              valuefile.printf(",\n  \"%s\" : {\n    \"action\" : \"%s\"",
                               actionBase->getLabel().c_str(),
                               actionBase->getName().c_str() );
            }
            if( avv ) {
              valuefile.printf(",\n    \"%s\" : { \"type\": \"atoms\","
                               " \"description\": \"virtual atom calculated by %s action\" }",
                               avv->getLabel().c_str(),
                               avv->getName().c_str() );
            } else {
              valuefile.printf(",\n    \"%s\" : { \"type\": \"atoms\","
                               " \"description\": \"indices of atoms specified in GROUP\" }",
                               actionBase->getLabel().c_str() );
            }
            valuefile.printf("\n  }");
            continue;
          }
          ActionWithValue* av=dynamic_cast<ActionWithValue*>(pp.get());
          if( av && av->getNumberOfComponents()>0 ) {
            Keywords keys;
            p.getKeywordsForAction( av->getName(), keys );
            if( firsta ) {
              valuefile.printf("  \"%s\" : {\n    \"action\" : \"%s\"",
                               av->getLabel().c_str(),
                               keys.getDisplayName().c_str() );
              firsta=false;
            } else {
              valuefile.printf(",\n  \"%s\" : {\n    \"action\" : \"%s\"",
                               av->getLabel().c_str(),
                               keys.getDisplayName().c_str() );
            }
            for(unsigned i=0; i<av->getNumberOfComponents(); ++i) {
              Value* myval = av->copyOutput(i);
              std::string compname = myval->getName(), description;
              if( av->getLabel()==compname ) {
                description = keys.getOutputComponentDescription(".#!value");
              } else {
                std::size_t dot=compname.find(av->getLabel() + ".");
                std::string cname = compname.substr(dot + av->getLabel().length() + 1);
                description = av->getOutputComponentDescription( cname, keys );
              }
              if( description.find("\\")!=std::string::npos ) {
                error("found invalid backslash character in documentation for component " + compname + " in action " + av->getName() + " with label " + av->getLabel() );
              }
              valuefile.printf(",\n    \"%s\" : { \"type\": \"%s\", \"description\": \"%s\" }", myval->getName().c_str(), myval->getValueType().c_str(), description.c_str() );
            }
            valuefile.printf("\n  }");
          }
          ActionShortcut* as=pp->castToActionShortcut();
          if( as ) {
            const auto & cnames = as->getSavedOutputs();
            if( cnames.size()==0 ) {
              continue ;
            }

            if( firsta ) {
              valuefile.printf("  \"shortcut_%s\" : {\n    \"action\" : \"%s\"", as->getShortcutLabel().c_str(), as->getName().c_str() );
              firsta=false;
            } else {
              valuefile.printf(",\n  \"shortcut_%s\" : {\n    \"action\" : \"%s\"", as->getShortcutLabel().c_str(), as->getName().c_str() );
            }
            Keywords keys;
            p.getKeywordsForAction( as->getName(), keys );
            for(unsigned i=0; i<cnames.size(); ++i) {
              ActionWithValue* av2=p.getActionSet().selectWithLabel<ActionWithValue*>( cnames[i] );
              if( !av2 ) {
                plumed_merror("could not find value created by shortcut with name " + cnames[i] );
              }
              if( av2->getNumberOfComponents()==1 ) {
                Value* myval = av2->copyOutput(0);
                std::string compname = myval->getName(), description;
                if( compname==as->getShortcutLabel() ) {
                  description = keys.getOutputComponentDescription(".#!value");
                } else {
                  std::size_t shortcutPos=compname.find(as->getShortcutLabel());
                  description = keys.getOutputComponentDescription(
                                  compname.substr(shortcutPos+as->getShortcutLabel().length()+1) );
                }
                if( description.find("\\")!=std::string::npos ) {
                  error("found invalid backslash character in documentation for component " + compname + " in action " + as->getName() + " with label " + as->getLabel() );
                }
                valuefile.printf(",\n    \"%s\" : { \"type\": \"%s\", \"description\": \"%s\" }",
                                 myval->getName().c_str(),
                                 myval->getValueType().c_str(),
                                 description.c_str() );
              } else {
                for(unsigned j=0; j<av2->getNumberOfComponents(); ++j) {
                  Value* myval = av2->copyOutput(j);
                  std::string compname = myval->getName(), description;
                  if( av2->getLabel()==compname ) {
                    plumed_merror("should not be outputting description of value from action when using shortcuts");
                  } else {
                    std::size_t dot=compname.find(av2->getLabel() + ".");
                    std::string cname = compname.substr(dot+av2->getLabel().length() + 1);
                    description = av2->getOutputComponentDescription( cname, keys );
                  }
                  if( description.find("\\")!=std::string::npos ) {
                    error("found invalid backslash character in documentation for component " + compname + " in action " + av2->getName() + " with label " + av2->getLabel() );
                  }
                  valuefile.printf(",\n    \"%s\" : { \"type\": \"%s\", \"description\": \"%s\" }", myval->getName().c_str(), myval->getValueType().c_str(), description.c_str() );
                }
              }
            }
            valuefile.printf("\n  }");
          }
        }
        valuefile.printf("\n}\n");
        valuefile.close();
      }
      if(parseOnly) {
        break;
      }
    }
    if(checknatoms!=natoms) {
      std::string stepstr;
      Tools::convert(step,stepstr);
      error("number of atoms in frame " + stepstr + " does not match number of atoms in first frame");
    }

    coordinates.assign(3*natoms,real(0.0));
    forces.assign(3*natoms,real(0.0));
    cell.assign(9,real(0.0));
    virial.assign(9,real(0.0));

    if( first_step || rnd.U01()>0.5) {
      if(debug_pd) {
        int npe=intracomm.Get_size();
        std::vector<int> loc(npe,0);
        std::vector<int> start(npe,0);
        for(int i=0; i<npe-1; i++) {
          int cc=(natoms*2*rnd.U01())/npe;
          if(start[i]+cc>natoms) {
            cc=natoms-start[i];
          }
          loc[i]=cc;
          start[i+1]=start[i]+loc[i];
        }
        loc[npe-1]=natoms-start[npe-1];
        intracomm.Bcast(loc,0);
        intracomm.Bcast(start,0);
        pd_nlocal=loc[intracomm.Get_rank()];
        pd_start=start[intracomm.Get_rank()];
        if(intracomm.Get_rank()==0) {
          std::fprintf(out,"\nDRIVER: Reassigning particle decomposition\n");
          std::fprintf(out,"DRIVER: ");
          for(int i=0; i<npe; i++) {
            std::fprintf(out,"%d ",loc[i]);
          }
          printf("\n");
          std::fprintf(out,"DRIVER: ");
          for(int i=0; i<npe; i++) {
            std::fprintf(out,"%d ",start[i]);
          }
          printf("\n");
        }
        p.cmd("setAtomsNlocal",pd_nlocal);
        p.cmd("setAtomsContiguous",pd_start);
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
          int n;
          for(n=0; n<npe; n++)
            if(n+1>r) {
              break;
            }
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
          std::fprintf(out,"\nDRIVER: Reassigning domain decomposition\n");
        }
        p.cmd("setAtomsNlocal",dd_nlocal);
        p.cmd("setAtomsGatindex",&dd_gatindex[0], {dd_nlocal});
      }
    }

    int plumedStopCondition=0;
    if(!noatoms) {
      auto errormessage=parser.readAtoms(
                          stride,
                          pbc_cli_given,
                          debug_pd,
                          pd_start,
                          pd_nlocal,
                          step,
                          masses.data(),
                          charges.data(),
                          coordinates.data(),
                          cell.data()
                        );

      if (errormessage) {
        if (*errormessage =="EOF") {
          break;
        } else {
          error(*errormessage);
        }
      }

      if(pbc_cli_given) {
        for(unsigned i=0; i<9; i++) {
          cell[i]=real(pbc_cli_box[i]);
        }
      }
      p.cmd("setStepLongLong",step);
      p.cmd("setStopFlag",&plumedStopCondition);

      if(debug_dd) {
        for(int i=0; i<dd_nlocal; ++i) {
          int kk=dd_gatindex[i];
          dd_coordinates[3*i+0]=coordinates[3*kk+0];
          dd_coordinates[3*i+1]=coordinates[3*kk+1];
          dd_coordinates[3*i+2]=coordinates[3*kk+2];
        }
        p.cmd("setForces",&dd_forces[0], {dd_nlocal,3});
        p.cmd("setPositions",&dd_coordinates[0], {dd_nlocal,3});
        p.cmd("setMasses",&dd_masses[0], {dd_nlocal});
        p.cmd("setCharges",&dd_charges[0], {dd_nlocal});
      } else {
// this is required to avoid troubles when the last domain
// contains zero atoms
// Basically, for empty domains we pass null pointers
#define fix_pd(xx) (pd_nlocal!=0?&xx:NULL)
        p.cmd("setForces",fix_pd(forces[3*pd_start]), {pd_nlocal,3});
        p.cmd("setPositions",fix_pd(coordinates[3*pd_start]), {pd_nlocal,3});
        p.cmd("setMasses",fix_pd(masses[pd_start]), {pd_nlocal});
        p.cmd("setCharges",fix_pd(charges[pd_start]), {pd_nlocal});
      }
      p.cmd("setBox",cell.data(), {3,3});
      p.cmd("setVirial",virial.data(), {3,3});
    } else {
      p.cmd("setStepLongLong",step);
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
    if(debug_pd) {
      intracomm.Sum(forces);
    }
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
        if(partner<0) {
          partner=0;
        }
        if(partner>=n) {
          partner=n-1;
        }
        p.cmd("GREX setPartner",partner);
        p.cmd("GREX calculate");
        p.cmd("GREX shareAllDeltaBias");
        for(int i=0; i<n; i++) {
          std::string s;
          Tools::convert(i,s);
          real a=std::numeric_limits<real>::quiet_NaN();
          s="GREX getDeltaBias "+s;
          p.cmd(s,&a);
          if(grex_log) {
            std::fprintf(grex_log," %f",a);
          }
        }
        if(grex_log) {
          std::fprintf(grex_log,"\n");
        }
      }
    }


    if(fp_forces) {
      std::fprintf(fp_forces,"%d\n",natoms);
      std::string fmtv=dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+"\n";
      std::string fmt=dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+"\n";
      if(dumpfullvirial) {
        std::fprintf(fp_forces,fmtv.c_str(),virial[0],virial[1],virial[2],virial[3],virial[4],virial[5],virial[6],virial[7],virial[8]);
      } else {
        std::fprintf(fp_forces,fmt.c_str(),virial[0],virial[4],virial[8]);
      }
      fmt="X "+fmt;
      for(int i=0; i<natoms; i++) {
        std::fprintf(fp_forces,fmt.c_str(),forces[3*i],forces[3*i+1],forces[3*i+2]);
      }
    }
    if(debugforces.length()>0) {
      // Now call the routine to work out the derivatives numerically
      numder.assign(3*natoms+9,real(0.0));
      real base=0;
      p.cmd("getBias",&base);
      if( std::abs(base)<epsilon ) {
        printf("WARNING: bias for configuration appears to be zero so debugging forces is trivial");
      }
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

    if(plumedStopCondition) {
      break;
    }

    step+=stride;
  }
  if(!parseOnly) {
    p.cmd("runFinalJobs");
  }

  return 0;
}

template<typename real>
void Driver<real>::evaluateNumericalDerivatives( const long long int& step, PlumedMain& p, const std::vector<real>& coordinates,
    const std::vector<real>& masses, const std::vector<real>& charges,
    std::vector<real>& cell, const double& base, std::vector<real>& numder ) {

  int natoms = coordinates.size() / 3;
  real delta = std::sqrt(epsilon);
  std::vector<Vector> pos(natoms);
  real bias=0;
  std::vector<real> fake_forces( 3*natoms ), fake_virial(9);
  for(int i=0; i<natoms; ++i) {
    for(unsigned j=0; j<3; ++j) {
      pos[i][j]=coordinates[3*i+j];
    }
  }

  for(int i=0; i<natoms; ++i) {
    for(unsigned j=0; j<3; ++j) {
      pos[i][j]=pos[i][j]+delta;
      p.cmd("setStepLongLong",step);
      p.cmd("setPositions",&pos[0][0], {natoms,3});
      p.cmd("setForces",&fake_forces[0], {natoms,3});
      p.cmd("setMasses",&masses[0], {natoms});
      p.cmd("setCharges",&charges[0], {natoms});
      p.cmd("setBox",&cell[0], {3,3});
      p.cmd("setVirial",&fake_virial[0], {3,3});
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
  Pbc pbc;
  pbc.setBox( box );
  Tensor nvirial;
  for(unsigned i=0; i<3; i++)
    for(unsigned k=0; k<3; k++) {
      double arg0=box(i,k);
      for(int j=0; j<natoms; ++j) {
        pos[j]=pbc.realToScaled( pos[j] );
      }
      cell[3*i+k]=box(i,k)=box(i,k)+delta;
      pbc.setBox(box);
      for(int j=0; j<natoms; j++) {
        pos[j]=pbc.scaledToReal( pos[j] );
      }
      p.cmd("setStepLongLong",step);
      p.cmd("setPositions",&pos[0][0], {natoms,3});
      p.cmd("setForces",&fake_forces[0], {natoms,3});
      p.cmd("setMasses",&masses[0], {natoms});
      p.cmd("setCharges",&charges[0], {natoms});
      p.cmd("setBox",&cell[0], {3,3});
      p.cmd("setVirial",&fake_virial[0], {3,3});
      p.cmd("prepareCalc");
      p.cmd("performCalcNoUpdate");
      p.cmd("getBias",&bias);
      cell[3*i+k]=box(i,k)=arg0;
      pbc.setBox(box);
      for(int j=0; j<natoms; j++)
        for(unsigned n=0; n<3; ++n) {
          pos[j][n]=coordinates[3*j+n];
        }
      nvirial(i,k) = ( bias - base ) / delta;
    }
  nvirial=-matmul(box.transpose(),nvirial);
  for(unsigned i=0; i<3; i++)
    for(unsigned k=0; k<3; k++) {
      numder[3*natoms+3*i+k] = nvirial(i,k);
    }

}

}
}
