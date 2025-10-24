/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "core/CLTool.h"
#include "core/CLToolRegister.h"
#include "tools/Tools.h"
#include "tools/Pbc.h"
#include "tools/PDB.h"
#include "core/ActionWithValue.h"
#include "core/ActionToPutData.h"
#include "core/ActionSet.h"
#include "core/Value.h"
#include "core/PlumedMain.h"
#include <cstdio>
#include <string>
#include <vector>

namespace PLMD {
namespace mapping {

//+PLUMEDOC TOOLS pathtools
/*
pathtools can be used to construct paths from pdb data

The path CVs in PLUMED are curvilinear coordinates through a high dimensional vector space.
Enhanced sampling calculations are often run using the progress along the paths and the distance from the path as CVs
as this provides a convenient way of defining a reaction coordinate for a complicated process.  This method is explained
in the documentation for [PATH](PATH.md).

The path itself is an ordered set of equally-spaced, high-dimensional frames the way in which these frames
should be constructed will depend on the problem in hand.  In other words, you will need to understand the reaction
you wish to study in order to select a sensible set of frames to use in your path CV.  This tool provides two
methods that may be useful when it comes to constructing paths; namely:

- A tool that takes in an initial guess path in which the frames are not equally spaced.  This tool adjusts the positions
of the frames in order to make them equally spaced so that they can be used as the basis for a path CV.

- A tool that takes two frames as input and that allows you to return a linear path connecting these two frames.  The
output from this method may be useful as an initial guess path.  It is arguable that a linear path rather defeats the
purpose of the path CV method, however, as the whole purpose is to be able to define non-linear paths.

## Examples

The example below shows how you can take a set of unequally spaced frames from a pdb file named `in_path.pdb`
and use pathtools to make them equally spaced so that they can be used as the basis for a path CV.  The file
containing this final path is named final_path.pdb.

```plumed
plumed pathtools --path in_path.pdb --metric EUCLIDEAN --out final_path.pdb
```

In this case, each frame in the path is defined using a set of collective variable values.  An extract from the input pdb file
will look something like this:

````
EMARK t1=-4.3345  t2=3.4725
END
REMARK t1=-4.1940  t2=3.3420
END
REMARK t1=-4.0535  t2=3.2114
END
REMARK t1=-3.9130  t2=3.0809
END
````

The example below shows how can create an initial linear path connecting the two pdb frames in start.pdb and
end.pdb.  In this case the path output to path.pdb will consist of 6 frames: the initial and final frames that
were contained in start.pdb and end.pdb as well as four equally spaced frames along the vector connecting
start.pdb to end.pdb.

```plumed
plumed pathtools --start start.pdb --end end.pdb --nframes 4 --metric OPTIMAL --out path.pdb
```

The vectors connecting frames here are calculated using the [RMSD](RMSD.md) action with the DISPLACEMENT option.  The positions on the path
are thus positions of atoms so the input files are normally formatted pdb files.

Often the idea with path collective variables is to create a path connecting some initial state A to some final state B.  You would
in this case have representative configurations from your A and B states defined in the input files to pathtools
that we have called `start.pdb` and `end.pdb` in the example above.  Furthermore, it may be useful to have a few frames
before your start frame and after your end frame.  You can use path tools to create these extended paths as shown below.
In this case the final path would now consist of 8 frames.  Four of these frames would lie on the vector connecting state
A to state B, there would be one frame each at start.pdb and end.pdb as well as one frame just before start.pdb and one
frame just after end.pdb.  All these frames would be equally spaced.

```plumed
plumed pathtools --start start.pdb --end end.pdb --nframes 4 --metric OPTIMAL --out path.pdb --nframes-before-start 2 --nframes-after-end 2
```

Notice also that when you re-parameterize paths you must choose two frames to fix.  Generally you chose to fix the states
that are representative of your states A and B.  By default pathtools will fix the first and last frames.  You can, however,
change the states to fix by taking advantage of the fixed flag as shown below.

```plumed
plumed pathtools --path inpath.pdb --metric EUCLIDEAN --out outpath.pdb --fixed 2,12
```

*/
//+ENDPLUMEDOC

class PathTools :
  public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit PathTools(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc);
  void printLambda( const std::string& mtype, const std::vector<std::string>& argstr, const std::string& ofile );
  std::string description()const {
    return "print out a description of the keywords for an action in html";
  }
};

PLUMED_REGISTER_CLTOOL(PathTools,"pathtools")

void PathTools::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("atoms","--start","a pdb file that contains the structure for the initial frame of your path");
  keys.add("atoms","--end","a pdb file that contains the structure for the final frame of your path");
  keys.add("atoms-1","--path","a pdb file that contains an initial path in which the frames are not equally spaced");
  keys.add("optional","--arg","the arguments that should be read in from the pdb files");
  keys.add("compulsory","--fixed","0","the frames to fix when constructing the path using --path");
  keys.add("compulsory","--metric","the measure to use to calculate the distance between frames");
  keys.add("compulsory","--out","the name of the file on which to output your path");
  keys.add("compulsory","--arg-fmt","%f","the format to use for argument values in your frames");
  keys.add("compulsory","--tolerance","1E-4","the tolerance to use for the algorithm that is used to re-parameterize the path");
  keys.add("compulsory","--nframes-before-start","1","the number of frames to include in the path before the first frame");
  keys.add("compulsory","--nframes","1","the number of frames between the start and end frames in your path");
  keys.add("compulsory","--nframes-after-end","1","the number of frames to put after the last frame of your path");
}

PathTools::PathTools(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}

int PathTools::main(FILE* in, FILE*out,Communicator& pc) {
  // Create a PLUMED object
  PlumedMain plmd;
  int s=sizeof(double);
  plmd.cmd("setRealPrecision",&s);
  plmd.cmd("setMDEngine","pathtools");
  // plmd.cmd("setNoVirial");
  double tstep=1.0;
  plmd.cmd("setTimestep",&tstep);
  plmd.cmd("init");
  int step=1;
  plmd.cmd("setStep",&step);

  std::string mtype;
  parse("--metric",mtype);
  std::string ifilename;
  parse("--path",ifilename);
  std::string ofmt;
  parse("--arg-fmt",ofmt);
  std::string ofilename;
  parse("--out",ofilename);
  std::vector<std::string> argstr;
  parseVector("--arg",argstr);
  if( ifilename.length()>0 ) {
    fprintf(out,"Reparameterising path in file named %s so that all frames are equally spaced \n",ifilename.c_str() );
    FILE* fp=fopen(ifilename.c_str(),"r");
    PDB mypdb;
    bool do_read=mypdb.readFromFilepointer(fp,false,0.1);
    std::vector<double> alig( mypdb.getOccupancy() ), disp( mypdb.getBeta() );
    std::vector<AtomNumber> indices( mypdb.getAtomNumbers() );
    std::vector<unsigned> residue_indices( indices.size() );
    for(unsigned i=0; i<residue_indices.size(); ++i) {
      residue_indices[i] = mypdb.getResidueNumber( indices[i] );
    }
    // Count the number of frames in the input file
    unsigned nfram=1;
    while ( do_read ) {
      if( !mypdb.readFromFilepointer(fp,false,0.1) ) {
        break;
      }
      nfram++;
    }
    if( argstr.size()>0 ) {
      for(unsigned i=0; i<argstr.size(); ++i ) {
        std::string input = argstr[i] + ": PDB2CONSTANT NOARGS REFERENCE=" + ifilename + " ARG=" + argstr[i];
        const char* inpt = input.c_str();
        plmd.cmd("readInputLine", inpt);
      }
    } else {
      std::string input = "data: PDB2CONSTANT REFERENCE=" + ifilename;
      const char* inpt = input.c_str();
      plmd.cmd("readInputLine", inpt );
    }
    std::string reparam_str = "REPARAMETERIZE_PATH REFERENCE=";
    if( argstr.size()>0 ) {
      reparam_str += argstr[0];
      for(unsigned i=0; i<argstr.size(); ++i ) {
        reparam_str += "," + argstr[i];
      }
    } else {
      reparam_str += "data";
    }
    std::vector<unsigned> fixed;
    parseVector("--fixed",fixed);
    if( fixed.size()==1 ) {
      if( fixed[0]!=0 ) {
        error("input to --fixed should be two integers");
      }
      fixed.resize(2);
      fixed[0]=1;
      fixed[1]=nfram;
    } else if( fixed.size()==2 ) {
      if( fixed[0]<1 || fixed[1]<1 || fixed[0]>nfram || fixed[1]>nfram ) {
        error("input to --fixed should be two numbers between 0 and the number of frames-1");
      }
    } else {
      error("input to --fixed should be two integers");
    }
    std::string fix1, fix2;
    Tools::convert( fixed[0], fix1 );
    Tools::convert( fixed[1], fix2 );
    reparam_str += " FIXED=" + fix1 + "," + fix2;
    std::string tol;
    parse("--tolerance",tol);
    reparam_str += " TOL=" + tol;
// Now create the metric object
    reparam_str += " METRIC=";
    if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) {
      reparam_str += "{RMSD DISPLACEMENT SQUARED UNORMALIZED TYPE=" + mtype;
      // Get the align values
      std::string anum;
      Tools::convert( alig[0], anum );
      reparam_str += " ALIGN=" + anum;
      for(unsigned i=1; i<alig.size(); ++i) {
        Tools::convert( alig[i], anum );
        reparam_str += "," + anum;
      }
      // Get the displace values
      std::string dnum;
      Tools::convert( disp[0], dnum );
      reparam_str += " DISPLACE=" + dnum;
      for(unsigned i=1; i<disp.size(); ++i) {
        Tools::convert( disp[i], dnum );
        reparam_str += "," + dnum;
      }
      reparam_str += "} METRIC_COMPONENT=disp";
    } else if( mtype=="EUCLIDEAN" ) {
      reparam_str += "DIFFERENCE";
    } else {
      // Add functionality to read plumed input here
      plumed_merror("metric type " + mtype + " has not been implemented");
    }

    // Now do the reparameterization
    const char* icinp= reparam_str.c_str();
    plmd.cmd("readInputLine",icinp);
    Action* raction = plmd.getActionSet()[plmd.getActionSet().size()-1].get();
    raction->update();

    // And print the final reference configurations
    std::string pinput="DUMPPDB STRIDE=1 DESCRIPTION=PATH FILE=" + ofilename + " FMT=" + ofmt;
    if( argstr.size()>0 ) {
      pinput += " ARG=" + argstr[0];
      for(unsigned i=1; i<argstr.size(); ++i ) {
        pinput += "," + argstr[i];
      }
    } else {
      std::string num;
      Tools::convert( indices[0].serial(), num );
      pinput += " ATOMS=data ATOM_INDICES=" + num;
      for(unsigned i=1; i<indices.size(); ++i ) {
        Tools::convert( indices[i].serial(), num );
        pinput += "," + num;
      }
      Tools::convert( residue_indices[0], num );
      pinput += " RESIDUE_INDICES=" + num;
      for(unsigned i=1; i<residue_indices.size(); ++i ) {
        Tools::convert( residue_indices[i], num );
        pinput += "," + num;
      }
      std::string anum, dnum;
      Tools::convert( alig[0], anum );
      Tools::convert( disp[0], dnum );
      pinput += " OCCUPANCY=" + anum;
      for(unsigned i=1; i<alig.size(); ++i) {
        Tools::convert( alig[i], anum );
        pinput += "," + anum;
      }
      pinput += " BETA=" + dnum;
      for(unsigned i=1; i<disp.size(); ++i) {
        Tools::convert( disp[i], dnum );
        pinput += "," + dnum;
      }
    }
    const char* pcinp=pinput.c_str();
    plmd.cmd("readInputLine",pcinp);
    Action* paction = plmd.getActionSet()[plmd.getActionSet().size()-1].get();
    paction->update();

    // Ouput data on spacings
    printLambda( mtype, argstr, ofilename );
    return 0;
  }
  for(unsigned i=0; i<argstr.size(); ++i) {
    std::string input = argstr[i] + ": CONSTANT VALUE=1";
    const char* inpt = input.c_str();
    plmd.cmd("readInputLine", inpt);
  }

// Read in the instructions
  unsigned nbefore, nbetween, nafter;
  std::string istart;
  parse("--start",istart);
  std::string iend;
  parse("--end",iend);
  parse("--nframes-before-start",nbefore);
  parse("--nframes",nbetween);
  parse("--nframes-after-end",nafter);
  nbetween++;
  fprintf(out,"Generating linear path connecting structure in file named %s to structure in file named %s \n",istart.c_str(),iend.c_str() );
  fprintf(out,"A path consisting of %u equally-spaced frames before the initial structure, %u frames between the initial and final structures "
          "and %u frames after the final structure will be created \n",nbefore,nbetween,nafter);

  std::vector<double> start_pos( argstr.size() ), end_pos( argstr.size() );
  std::vector<AtomNumber> indices;
  if( argstr.size()>0 ) {
    for(unsigned i=0; i<argstr.size(); ++i ) {
      std::string input = argstr[i] + "_start: PDB2CONSTANT REFERENCE=" + istart + " ARG=" + argstr[i];
      const char* inpt = input.c_str();
      plmd.cmd("readInputLine", inpt );
      long srank;
      plmd.cmd("getDataRank " + argstr[i] + "_start", &srank );
      if( srank!=1 ) {
        error("should only be one frame in input pdb");
      }
      std::vector<long> sshape( 1 );
      plmd.cmd("getDataShape " + argstr[i] + "_start", &sshape[0] );
      if( sshape[0]!=1 ) {
        error("should only be one frame in input pdb");
      }
      plmd.cmd("setMemoryForData " + argstr[i] + "_start", &start_pos[i] );
      std::string input2 = argstr[i] + "_end: PDB2CONSTANT REFERENCE=" + iend + " ARG=" + argstr[i];
      const char* inpt2 = input2.c_str();
      plmd.cmd("readInputLine", inpt2 );
      long erank;
      plmd.cmd("getDataRank " + argstr[i] + "_end", &erank );
      if( erank!=1 ) {
        error("should only be one frame in input pdb");
      }
      std::vector<long> eshape( 1 );
      plmd.cmd("getDataShape " + argstr[i] + "_end", &eshape[0] );
      if( eshape[0]!=1 ) {
        error("should only be one frame in input pdb");
      }
      plmd.cmd("setMemoryForData " + argstr[i] + "_end", &end_pos[i] );
    }
  } else {
    std::string input = "start: PDB2CONSTANT REFERENCE=" + istart;
    const char* inpt = input.c_str();
    plmd.cmd("readInputLine", inpt );
    FILE* fp=fopen(istart.c_str(),"r");
    PDB mypdb;
    mypdb.readFromFilepointer(fp,false,0.1);
    indices.resize( mypdb.getAtomNumbers().size() );
    for(unsigned i=0; i<indices.size(); ++i) {
      indices[i] = mypdb.getAtomNumbers()[i];
    }
    long srank;
    plmd.cmd("getDataRank start", &srank );
    if( srank!=2 ) {
      error("should only be one frame in input pdb:" + istart);
    }
    std::vector<long> sshape( srank );
    plmd.cmd("getDataShape start", &sshape[0] );
    start_pos.resize( sshape[0]*sshape[1] );
    if( sshape[0]!=1 ) {
      error("should only be one frame in input pdb: " + istart);
    }
    plmd.cmd("setMemoryForData start", &start_pos[0] );
    std::string input2 = "end: PDB2CONSTANT REFERENCE=" + iend;
    const char* inpt2 = input2.c_str();
    plmd.cmd("readInputLine", inpt2 );
    long erank;
    plmd.cmd("getDataRank end", &erank );
    if( erank!=2 ) {
      error("should only be one frame in input pdb: " + iend);
    }
    std::vector<long> eshape( erank );
    plmd.cmd("getDataShape end", &eshape[0] );
    end_pos.resize( eshape[0]*eshape[1] );
    if( eshape[0]!=1 ) {
      error("should only be one frame in input pdb: " + iend );
    }
    plmd.cmd("setMemoryForData end", &end_pos[0] );
  }

// Now create the metric object
  std::vector<double> alig, disp;
  std::vector<unsigned> residue_indices;
  if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) {
    std::string trinput = "endT: TRANSPOSE ARG=end";
    const char* tinp=trinput.c_str();
    plmd.cmd("readInputLine",tinp);
    PDB pdb;
    pdb.read(istart,false,0.1);
    alig.resize( pdb.getOccupancy().size() );
    disp.resize( pdb.getBeta().size() );
    for(unsigned i=0; i<alig.size(); ++i) {
      alig[i] = pdb.getOccupancy()[i];
      disp[i]= pdb.getBeta()[i];
    }
    std::string minput = "d: RMSD DISPLACEMENT SQUARED UNORMALIZED TYPE=" + mtype + " ARG=endT,start";
    // Get the align values
    std::string anum;
    Tools::convert( alig[0], anum );
    minput += " ALIGN=" + anum;
    for(unsigned i=1; i<alig.size(); ++i) {
      Tools::convert( alig[i], anum );
      minput += "," + anum;
    }
    // Get the displace values
    std::string dnum;
    Tools::convert( disp[0], dnum );
    minput += " DISPLACE=" + dnum;
    for(unsigned i=1; i<disp.size(); ++i) {
      Tools::convert( disp[i], dnum );
      minput += "," + dnum;
    }
    const char* mcinp=minput.c_str();
    plmd.cmd("readInputLine",mcinp);
    std::string tinput = "CUSTOM ARG=d.disp FUNC=x PERIODIC=NO";
    const char* tcinp=tinput.c_str();
    plmd.cmd("readInputLine",tcinp);
    // Get the residue numbers
    residue_indices.resize( pdb.size() );
    for(unsigned i=0; i<residue_indices.size(); ++i) {
      residue_indices[i] = pdb.getResidueNumber( indices[i] );
    }
  } else if( mtype=="EUCLIDEAN" ) {
    std::string end_args="ARG1=" + argstr[0] + "_end", start_args="ARG2=" + argstr[0] + "_start";
    for(unsigned i=1; i<argstr.size(); ++i ) {
      end_args += "," + argstr[i] + "_end";
      start_args += "," + argstr[i] + "_start";
    }
    std::string minput = "d: DISPLACEMENT " + end_args + " " + start_args;
    const char* mcinp=minput.c_str();
    plmd.cmd("readInputLine",mcinp);
    std::string tinput = "TRANSPOSE ARG=d";
    const char* tcinp=tinput.c_str();
    plmd.cmd("readInputLine",tcinp);
  } else {
    // Add functionality to read plumed input here
    plumed_merror("metric type " + mtype + " has not been implemented");
  }

// Retrieve final displacement vector
  unsigned aind = plmd.getActionSet().size()-1;
  while( true ) {
    const ActionShortcut* as=dynamic_cast<const ActionShortcut*>( plmd.getActionSet()[aind].get() );
    if( !as ) {
      break ;
    }
    aind = aind - 1;
    plumed_assert( aind>=0 );
  }
  ActionWithValue* av = dynamic_cast<ActionWithValue*>( plmd.getActionSet()[aind].get() );
  if( !av ) {
    error("invalid input for metric" );
  }
  if( av->getNumberOfComponents()!=1 && av->getName()!="RMSD" ) {
    error("cannot use multi component actions as metric");
  }
  std::string mydisp = av->copyOutput(0)->getName();
// Now add calls so we can grab the data from plumed
  long rank;
  plmd.cmd("getDataRank " + mydisp, &rank );
  if( rank!=1 ) {
    error("displacement must be a vector quantity");
  }
  std::vector<long> ishape( rank );
  plmd.cmd("getDataShape " + mydisp, &ishape[0] );
  std::vector<double> displacement( ishape[0] );
  plmd.cmd("setMemoryForData " + mydisp, &displacement[0] );
// And calculate the displacement
  plmd.cmd("calc");

  // Now read in the initial frame

  // Now create frames
  double delr = 1.0 / static_cast<double>( nbetween );
  std::vector<std::vector<double> > allframes;
  std::vector<double> frame( displacement.size() );
  for(unsigned i=0; i<nbefore; ++i) {
    for(unsigned j=0; j<frame.size(); ++j) {
      frame[j] = start_pos[j] - i*delr*displacement[j];
    }
    allframes.push_back( frame );
  }
  for(unsigned i=1; i<nbetween; ++i) {
    for(unsigned j=0; j<frame.size(); ++j) {
      frame[j] = start_pos[j] + i*delr*displacement[j];
    }
    allframes.push_back( frame );
  }
  for(unsigned i=0; i<nafter; ++i) {
    for(unsigned j=0; j<frame.size(); ++j) {
      frame[j] = end_pos[j] + i*delr*displacement[j];
    }
    allframes.push_back( frame );
  }
  // This prints out our final reference configurations
  plmd.cmd("clear");
  plmd.readInputLine("timestep: PUT UNIT=time PERIODIC=NO CONSTANT");
  ActionToPutData* ts = plmd.getActionSet().selectWithLabel<ActionToPutData*>("timestep");
  ts->setValuePointer("timestep", &tstep );
  std::string pinput="DUMPPDB STRIDE=1 DESCRIPTION=PATH FILE=" + ofilename + " FMT=" + ofmt;
  if( argstr.size()>0 ) {
    for(unsigned i=0; i<argstr.size(); ++i) {
      std::string rnum;
      Tools::convert( allframes[0][i], rnum );
      std::string inpt = argstr[i] + ": CONSTANT VALUES=" + rnum;
      for(unsigned j=1; j<allframes.size(); ++j) {
        Tools::convert( allframes[j][i], rnum );
        inpt += "," + rnum;
      }
      const char* icinp=inpt.c_str();
      plmd.cmd("readInputLine",icinp);
      if( i==0 ) {
        pinput += " ARG=" + argstr[i];
      } else {
        pinput += "," + argstr[i];
      }
    }
  } else {
    std::string nc;
    Tools::convert( frame.size(), nc );
    std::string nr;
    Tools::convert( allframes.size(), nr );
    std::string rnum;
    Tools::convert( allframes[0][0], rnum );
    std::string inpt = "atom_data: CONSTANT NROWS=" + nr + " NCOLS=" + nc + " VALUES=" + rnum;
    for(unsigned i=0; i<allframes.size(); ++i) {
      for(unsigned j=0; j<allframes[i].size(); ++j) {
        if( i==0 && j==0 ) {
          continue;
        }
        Tools::convert( allframes[i][j], rnum );
        inpt += "," + rnum;
      }
    }
    const char* icinp=inpt.c_str();
    plmd.cmd("readInputLine",icinp);
    std::string num;
    Tools::convert( indices[0].serial(), num );
    pinput += " ATOMS=atom_data ATOM_INDICES=" + num;
    for(unsigned i=1; i<indices.size(); ++i ) {
      Tools::convert( indices[i].serial(), num );
      pinput += "," + num;
    }
    Tools::convert( residue_indices[0], num );
    pinput += " RESIDUE_INDICES=" + num;
    for(unsigned i=1; i<residue_indices.size(); ++i ) {
      Tools::convert( residue_indices[i], num );
      pinput += "," + num;
    }
    std::string anum, dnum;
    Tools::convert( alig[0], anum );
    Tools::convert( disp[0], dnum );
    pinput += " OCCUPANCY=" + anum;
    for(unsigned i=1; i<alig.size(); ++i) {
      Tools::convert( alig[i], anum );
      pinput += "," + anum;
    }
    pinput += " BETA=" + dnum;
    for(unsigned i=1; i<disp.size(); ++i) {
      Tools::convert( disp[i], dnum );
      pinput += "," + dnum;
    }
  }
  const char* pcinp=pinput.c_str();
  plmd.cmd("readInputLine",pcinp);
  Action* paction = plmd.getActionSet()[plmd.getActionSet().size()-1].get();
  paction->update();
  // And output suggestions on the value of Lambda
  printLambda( mtype, argstr, ofilename );
  return 0;
}

void PathTools::printLambda( const std::string& mtype, const std::vector<std::string>& argstr, const std::string& ofile ) {
  // Create a PLUMED object
  PlumedMain plmd;
  int s=sizeof(double);
  plmd.cmd("setRealPrecision",&s);
  plmd.cmd("setMDEngine","pathtools");
  double tstep=1.0;
  plmd.cmd("setTimestep",&tstep);
  plmd.cmd("init");
  int step=1;
  plmd.cmd("setStep",&step);

  FILE* fp=fopen(ofile.c_str(),"r");
  bool do_read=true;
  unsigned nfram=0;
  std::vector<double> alig, disp;
  // Read in the argument names
  for(unsigned i=0; i<argstr.size(); ++i ) {
    std::string input2 = argstr[i] + ": CONSTANT VALUE=1";
    const char* inpt2 = input2.c_str();
    plmd.cmd("readInputLine", inpt2);
  }
  while (do_read ) {
    PDB mypdb;
    // Read the pdb file
    do_read=mypdb.readFromFilepointer(fp,false,0.1);
    if( !do_read ) {
      break;
    }
    std::string num;
    Tools::convert( nfram+1, num );
    nfram++;
    std::string iinput;
    if( argstr.size()>0 ) {
      for(unsigned i=0; i<argstr.size(); ++i ) {
        std::string input = "ref_" + num + "_" + argstr[i] + ": PDB2CONSTANT REFERENCE=" + ofile + " ARG=" + argstr[i] + " NUMBER=" + num;
        const char* inpt = input.c_str();
        plmd.cmd("readInputLine", inpt );
      }
    } else {
      std::string input = "ref_" + num + ": PDB2CONSTANT REFERENCE=" + ofile + " NUMBER=" + num;
      const char* inpt = input.c_str();
      plmd.cmd("readInputLine", inpt );
    }

    if( nfram==1 ) {
      alig.resize( mypdb.getOccupancy().size() );
      for(unsigned i=0; i<alig.size(); ++i) {
        alig[i]=mypdb.getOccupancy()[i];
      }
      disp.resize( mypdb.getBeta().size() );
      for(unsigned i=0; i<disp.size(); ++i) {
        disp[i]=mypdb.getBeta()[i];
      }
    }
  }
  // Now create the objects to measure the distances between the frames
  std::vector<double> data( nfram );
  for(unsigned j=1; j<nfram; ++j) {
    std::string istr, jstr;
    Tools::convert( j, istr);
    Tools::convert( j+1, jstr );
    if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) {
      std::string strv = "ref_" + jstr + "T: TRANSPOSE ARG=ref_" + jstr;
      const char* cstrv = strv.c_str();
      plmd.cmd("readInputLine",cstrv);
      std::string dstring = "d" + istr + ": RMSD TYPE=" + mtype + " ARG=ref_" + jstr + "T,ref_" + istr;
      // Get the align values
      std::string anum;
      Tools::convert( alig[0], anum );
      dstring += " ALIGN=" + anum;
      for(unsigned i=1; i<alig.size(); ++i) {
        Tools::convert( alig[i], anum );
        dstring += "," + anum;
      }
      // Get the displace values
      std::string dnum;
      Tools::convert( disp[0], dnum );
      dstring += " DISPLACE=" + dnum;
      for(unsigned i=1; i<disp.size(); ++i) {
        Tools::convert( disp[i], dnum );
        dstring += "," + dnum;
      }
      const char* icinp=dstring.c_str();
      plmd.cmd("readInputLine",icinp);
    } else if( mtype=="EUCLIDEAN" ) {
      std::string end_args="ARG1=ref_" + istr + "_" + argstr[0], start_args="ARG2=ref_" + jstr + "_" + argstr[0];
      for(unsigned i=1; i<argstr.size(); ++i ) {
        end_args += ",ref_" + istr + "_" + argstr[i];
        start_args += ",ref_" + jstr + "_" + argstr[i];
      }
      std::string fstr = "d" + istr + ": EUCLIDEAN_DISTANCE " + end_args + " " + start_args;
      const char* icinp=fstr.c_str();
      plmd.cmd("readInputLine",icinp);
    } else {
      plumed_merror("metric type " + mtype + " has not been implemented");
    }
    long rank;
    plmd.cmd("getDataRank d" + istr, &rank );
    if( rank!=1 ) {
      error("distance should be of rank 1");
    }
    std::vector<long> ishape(1);
    plmd.cmd("getDataShape d" + istr, &ishape[0] );
    if( ishape[0]!=1 ) {
      error("distance should be of rank 1");
    }
    plmd.cmd("setMemoryForData d" + istr, &data[j-1] );
  }
  plmd.cmd("calc");
  double mean=0;
  printf("N.B. THIS CODE ALWAYS AIMS TO CREATE EQUALLY SPACED FRAMES \n");
  printf("THERE MAY BE SMALL DESCREPENCIES IN THE NUMBERS BELOW, HOWEVER, BECAUSE OF ROUNDING ERRORS \n");
  for(unsigned i=1; i<nfram; ++i) {
    printf("FINAL DISTANCE BETWEEN FRAME %u AND %u IS %f \n",i-1,i,data[i-1]);
    mean += data[i-1];
  }
  printf("SUGGESTED LAMBDA PARAMETER IS THUS %f \n",2.3/mean/static_cast<double>( nfram-1 ) );
}

} // End of namespace
}
