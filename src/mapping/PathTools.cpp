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
#include "cltools/CLTool.h"
#include "cltools/CLToolRegister.h"
#include "tools/Tools.h"
#include "tools/Pbc.h"
#include "tools/PDB.h"
#include "core/ActionWithValue.h"
#include "core/ActionSet.h"
#include "core/Value.h"
#include "core/PlumedMain.h"
#include "setup/SetupReferenceBase.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

namespace PLMD {
namespace mapping {

//+PLUMEDOC TOOLS pathtools
/*
pathtools can be used to construct paths from pdb data

The path CVs in PLUMED are curvilinear coordinates through a high dimensional vector space.
Enhanced sampling calculations are often run using the progress along the paths and the distance from the path as CVs
as this provides a convenient way of defining a reaction coordinate for a complicated process.  This method is explained
in the documentation for \ref PATH.

The path itself is an ordered set of equally-spaced, high-dimensional frames the way in which these frames
should be constructed will depend on the problem in hand.  In other words, you will need to understand the reaction
you wish to study in order to select a sensible set of frames to use in your path CV.  This tool provides two
methods that may be useful when it comes to constructing paths; namely:

- A tool that takes in an initial guess path in which the frames are not equally spaced.  This tool adjusts the positions
of the frames in order to make them equally spaced so that they can be used as the basis for a path CV.

- A tool that takes two frames as input and that allows you to return a linear path connecting these two frames.  The
output from this method may be useful as an initial guess path.  It is arguable that a linear path rather defeats the
purpose of the path CV method, however, as the whole purpose is to be able to define non-linear paths.

Notice that you can use these two methods and take advantage of all the ways of measuring \ref dists that are available within
PLUMED. The way you do this with each of these tools described above is explained in the example below.

\par Examples

The example below shows how you can take a set of unequally spaced frames from a pdb file named in_path.pdb
and use pathtools to make them equally spaced so that they can be used as the basis for a path CV.  The file
containing this final path is named final_path.pdb.

\verbatim
plumed pathtools --path in_path.pdb --metric EUCLIDEAN --out final_path.pdb
\endverbatim

The example below shows how can create an initial linear path connecting the two pdb frames in start.pdb and
end.pdb.  In this case the path output to path.pdb will consist of 6 frames: the initial and final frames that
were contained in start.pdb and end.pdb as well as four equally spaced frames along the vector connecting
start.pdb to end.pdb.

\verbatim
plumed pathtools --start start.pdb --end end.pdb --nframes 4 --metric OPTIMAL --out path.pdb
\endverbatim

Often the idea with path collective variables is to create a path connecting some initial state A to some final state B.  You would
in this case have representative configurations from your A and B states defined in the input files to pathtools
that we have called start.pdb and end.pdb in the example above.  Furthermore, it may be useful to have a few frames
before your start frame and after your end frame.  You can use path tools to create these extended paths as shown below.
In this case the final path would now consist of 8 frames.  Four of these frames would lie on the vector connecting state
A to state B, there would be one frame each at start.pdb and end.pdb as well as one frame just before start.pdb and one
frame just after end.pdb.  All these frames would be equally spaced.

\verbatim
plumed pathtools --start start.pdb --end end.pdb --nframes 4 --metric OPTIMAL --out path.pdb --nframes-before-start 2 --nframes-after-end 2
\endverbatim

Notice also that when you re-parameterize paths you must choose two frames to fix.  Generally you chose to fix the states
that are representative of your states A and B.  By default pathtools will fix the first and last frames.  You can, however,
change the states to fix by taking advantage of the fixed flag as shown below.

\verbatim
plumed pathtools --path inpath.pdb --metric EUCLIDEAN --out outpath.pdb --fixed 2,12
\endverbatim

*/
//+ENDPLUMEDOC

class PathTools :
  public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  explicit PathTools(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc);
  void printLambda( const std::string& mtype, const std::string& argstr, const std::string& ofile );
  string description()const {
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
  CLTool(co)
{
  inputdata=commandline;
}

int PathTools::main(FILE* in, FILE*out,Communicator& pc) {
  // Create a PLUMED object
  PlumedMain plmd; int s=sizeof(double);
  plmd.cmd("setRealPrecision",&s);
  plmd.cmd("setNoVirial"); 
  plmd.cmd("setMDEngine","pathtools");
  int natoms = 1; plmd.cmd("setNatoms",&natoms);
  double tstep=1.0; plmd.cmd("setTimestep",&tstep);
  plmd.cmd("init");
  int step=1; plmd.cmd("setStep",&step);
  std::vector<double> mass(1); plmd.cmd("setMasses",&mass[0]); 
  std::vector<double> charge(1); plmd.cmd("setCharges",&charge[0]);
  std::vector<Vector> pos(1); plmd.cmd("setPositions",&pos[0][0]); 
  std::vector<Vector> forces(1); plmd.cmd("setForces",&forces[0][0]);

  std::string mtype; parse("--metric",mtype);
  std::string ifilename; parse("--path",ifilename);
  std::string ofmt; parse("--arg-fmt",ofmt);
  std::string ofilename; parse("--out",ofilename);
  std::string argstr; parse("--arg",argstr);
  if( ifilename.length()>0 ) {
    fprintf(out,"Reparameterising path in file named %s so that all frames are equally spaced \n",ifilename.c_str() );
    FILE* fp=fopen(ifilename.c_str(),"r");
    bool do_read=true; unsigned nfram=0;
    std::vector<double> alig, disp; std::vector<AtomNumber> indices;
    while (do_read) {
      PDB mypdb;
      // Read the pdb file
      do_read=mypdb.readFromFilepointer(fp,false,0.1);
      if( do_read ) {
          std::string num; Tools::convert( nfram+1, num ); nfram++;
          std::string iinput = "ref_" + num + ": READ_CONFIG REFERENCE=" + ifilename + " NUMBER=" + num;
          if( argstr.length()>0 ) iinput += " READ_ARG=" + argstr;
          const char* icinp=iinput.c_str(); plmd.cmd("readInputLine",icinp);
          if( nfram==1 ) {
              indices.resize( mypdb.getAtomNumbers().size() );
              for(unsigned i=0;i<indices.size();++i) indices[i]=mypdb.getAtomNumbers()[i];
              alig.resize( mypdb.getOccupancy().size() );
              for(unsigned i=0;i<alig.size();++i) alig[i]=mypdb.getOccupancy()[i];
              disp.resize( mypdb.getBeta().size() );
              for(unsigned i=0;i<disp.size();++i) disp[i]=mypdb.getBeta()[i];
          }
      }
    }
    std::string reparam_str = "REPARAMETERIZE_PATH REFFRAMES=ref_1"; 
    for(unsigned i=1;i<nfram;++i){ std::string num; Tools::convert(i+1,num); reparam_str += ",ref_" + num; } 
    std::vector<unsigned> fixed; parseVector("--fixed",fixed);
    if( fixed.size()==1 ) {
      if( fixed[0]!=0 ) error("input to --fixed should be two integers");
      fixed.resize(2); fixed[0]=1; fixed[1]=nfram;
    } else if( fixed.size()==2 ) {
      if( fixed[0]<1 || fixed[1]<1 || fixed[0]>nfram || fixed[1]>nfram ) {
        error("input to --fixed should be two numbers between 0 and the number of frames-1");
      }
    } else error("input to --fixed should be two integers");
    std::string fix1, fix2; Tools::convert( fixed[0], fix1 ); Tools::convert( fixed[1], fix2 );
    reparam_str += " FIXED=" + fix1 + "," + fix2;
    std::string tol; parse("--tolerance",tol); reparam_str += " TOL=" + tol; 
// Now create the metric object
    reparam_str += " METRIC={";
    if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) {
        std::string atnum; Tools::convert( indices[0].serial(), atnum ); 
        reparam_str += "RMSD DISPLACEMENT SQUARED UNORMALIZED TYPE=" + mtype + " REFERENCE_ATOMS=" + atnum;
        for(unsigned i=1;i<indices.size();++i){ Tools::convert( indices[i].serial(), atnum ); reparam_str += "," + atnum; }       
        natoms=indices[0].serial(); 
        for(unsigned i=1;i<indices.size();++i) {
          if( indices[i].serial()>natoms ) natoms = indices[i].serial();
        }
        Tools::convert( natoms+indices[0].serial(), atnum ); reparam_str += " ATOMS=" + atnum;
        for(unsigned i=1;i<alig.size();++i){ Tools::convert( natoms + indices[i].serial(), atnum ); reparam_str += "," + atnum; }  
        // Get the align values 
        std::string anum; Tools::convert( alig[0], anum ); reparam_str += " ALIGN=" + anum;
        for(unsigned i=1;i<alig.size();++i){ Tools::convert( alig[i], anum ); reparam_str += "," + anum; }
        // Get the displace values
        std::string dnum; Tools::convert( disp[0], dnum ); reparam_str += " DISPLACE=" + dnum;
        for(unsigned i=1;i<disp.size();++i){ Tools::convert( disp[i], dnum ); reparam_str += "," + dnum; }
    } else if( mtype=="EUCLIDEAN" ) {
        reparam_str += "DIFFERENCE ARG1=arg2 ARG2=arg1";
    } else {
       // Add functionality to read plumed input here
       plumed_merror("metric type " + mtype + " has not been implemented");
    }
    reparam_str += "}";

    // Now do the reparameterization
    const char* icinp= reparam_str.c_str(); plmd.cmd("readInputLine",icinp);
    Action* raction = plmd.getActionSet()[plmd.getActionSet().size()-1].get();
    raction->update();

    // And print the final reference configurations
    std::string pinput="PRINT STRIDE=1 DESCRIPTION=PATH FILE=" + ofilename + " FMT=" + ofmt;
    for(unsigned i=0;i<nfram;++i){ std::string num; Tools::convert( i+1, num ); pinput += " CONFIG" + num + "=ref_" + num; }
    const char* pcinp=pinput.c_str(); plmd.cmd("readInputLine",pcinp);
    Action* paction = plmd.getActionSet()[plmd.getActionSet().size()-1].get();
    ActionAtomistic* aact = dynamic_cast<ActionAtomistic*>( paction );
    aact->retrieveAtoms(); paction->update();

    // Ouput data on spacings
    printLambda( mtype, argstr, ofilename );
    return 0;
  }

// Read in the instructions
  unsigned nbefore, nbetween, nafter;
  std::string istart; parse("--start",istart); std::string iend; parse("--end",iend); 
  parse("--nframes-before-start",nbefore); parse("--nframes",nbetween); parse("--nframes-after-end",nafter);
  nbetween++;
  fprintf(out,"Generating linear path connecting structure in file named %s to structure in file named %s \n",istart.c_str(),iend.c_str() );
  fprintf(out,"A path consisting of %u equally-spaced frames before the initial structure, %u frames between the intial and final structures "
          "and %u frames after the final structure will be created \n",nbefore,nbetween,nafter);

// Read initial frame
  std::string iinput = "start: READ_CONFIG REFERENCE=" + istart;
  if( argstr.length()>0 ) iinput += " READ_ARG=" + argstr;
  const char* icinp=iinput.c_str(); plmd.cmd("readInputLine",icinp);

// Read final frame
  std::string finput = "end: READ_CONFIG REFERENCE=" + iend;
  if( argstr.length()>0 ) finput += " READ_ARG=" + argstr;
  const char* fcinp=finput.c_str(); plmd.cmd("readInputLine",fcinp);

// Now create the metric object
  if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) { 
      PDB pdb; pdb.read(istart,false,0.1); 
      std::vector<double> alig( pdb.getOccupancy() ), disp( pdb.getBeta() );
      std::string minput = "RMSD DISPLACEMENT SQUARED UNORMALIZED TYPE=" + mtype + " REFERENCE_ATOMS=start ATOMS=end";
      // Get the align values 
      std::string anum; Tools::convert( alig[0], anum ); minput += " ALIGN=" + anum;
      for(unsigned i=1;i<alig.size();++i){ Tools::convert( alig[i], anum ); minput += "," + anum; }
      // Get the displace values
      std::string dnum; Tools::convert( disp[0], dnum ); minput += " DISPLACE=" + dnum;
      for(unsigned i=1;i<disp.size();++i){ Tools::convert( disp[i], dnum ); minput += "," + dnum; }   
      const char* mcinp=minput.c_str(); plmd.cmd("readInputLine",mcinp); 
  } else if( mtype=="EUCLIDEAN" ) {
      std::string minput = "DIFFERENCE ARG1=end ARG2=start";
      const char* mcinp=minput.c_str(); plmd.cmd("readInputLine",mcinp);
  } else {
     // Add functionality to read plumed input here
     plumed_merror("metric type " + mtype + " has not been implemented");
  }

// Retrieve final displacement vector
  ActionWithValue* av=dynamic_cast<ActionWithValue*>( plmd.getActionSet()[plmd.getActionSet().size()-1].get() );
  if( !av ) error("invalid input for metric" );
  if( av->getNumberOfComponents()!=1 && av->getName()!="RMSD" ) error("cannot use multi component actions as metric");
  std::string mydisp = av->copyOutput(0)->getName();
// Now add calls so we can grab the data from plumed
  long rank; plmd.cmd("getDataRank " + mydisp, &rank );
  if( rank!=1 ) error("displacement must be a vector quantity");
  std::vector<long> ishape( rank ); plmd.cmd("getDataShape " + mydisp, &ishape[0] );
  std::vector<double> displacement( ishape[0] ); plmd.cmd("setMemoryForData " + mydisp, &displacement[0] );
// And calculate the displacement
  plmd.cmd("calc");

  // Now create frames
  double delr = 1.0 / static_cast<double>( nbetween ); unsigned nframes=0;
  for(int i=0; i<nbefore; ++i) {
      std::string num; Tools::convert( nframes+1, num ); nframes++;
      std::string fstr = "frame" + num + ": READ_CONFIG REFERENCE=" + istart;
      if( argstr.length()>0 ) fstr += " READ_ARG=" + argstr;
      const char* ifstr=fstr.c_str(); plmd.cmd("readInputLine",ifstr); 
      setup::SetupReferenceBase* sb = dynamic_cast<setup::SetupReferenceBase*>( plmd.getActionSet()[plmd.getActionSet().size()-1].get() );
      sb->displaceReferenceConfiguration( -i*delr, displacement );
  }
  for(unsigned i=1; i<nbetween; ++i) {
    std::string num; Tools::convert( nframes+1, num ); nframes++;
    std::string fstr = "frame" + num + ": READ_CONFIG REFERENCE=" + istart;
    if( argstr.length()>0 ) fstr += " READ_ARG=" + argstr;
    const char* ifstr=fstr.c_str(); plmd.cmd("readInputLine",ifstr); 
    setup::SetupReferenceBase* sb = dynamic_cast<setup::SetupReferenceBase*>( plmd.getActionSet()[plmd.getActionSet().size()-1].get() );
    sb->displaceReferenceConfiguration( i*delr, displacement );
  }
  for(unsigned i=0; i<nafter; ++i) {
    std::string num; Tools::convert( nframes+1, num ); nframes++;
    std::string fstr = "frame" + num + ": READ_CONFIG REFERENCE=" + iend;
    if( argstr.length()>0 ) fstr += " READ_ARG=" + argstr;
    const char* ifstr=fstr.c_str(); plmd.cmd("readInputLine",ifstr);
    setup::SetupReferenceBase* sb = dynamic_cast<setup::SetupReferenceBase*>( plmd.getActionSet()[plmd.getActionSet().size()-1].get() );
    sb->displaceReferenceConfiguration( i*delr, displacement );
  }

  // This prints out our final reference configurations
  std::string pinput="PRINT STRIDE=1 DESCRIPTION=PATH FILE=" + ofilename + " FMT=" + ofmt; 
  for(unsigned i=0;i<nframes;++i){ std::string num; Tools::convert( i+1, num ); pinput += " CONFIG" + num + "=frame" + num; }
  const char* pcinp=pinput.c_str(); plmd.cmd("readInputLine",pcinp); 
  Action* paction = plmd.getActionSet()[plmd.getActionSet().size()-1].get();
  ActionAtomistic* aact = dynamic_cast<ActionAtomistic*>( paction ); 
  aact->retrieveAtoms(); paction->update();
  // And output suggestions on the value of Lambda
  printLambda( mtype, argstr, ofilename );
  return 0;
}

void PathTools::printLambda( const std::string& mtype, const std::string& argstr, const std::string& ofile ) {
  // Create a PLUMED object
  PlumedMain plmd; int s=sizeof(double);
  plmd.cmd("setRealPrecision",&s);
  plmd.cmd("setNoVirial");
  plmd.cmd("setMDEngine","pathtools");
  int natoms = 1; plmd.cmd("setNatoms",&natoms);
  double tstep=1.0; plmd.cmd("setTimestep",&tstep);
  plmd.cmd("init");
  int step=1; plmd.cmd("setStep",&step);
  std::vector<double> mass(1); plmd.cmd("setMasses",&mass[0]); 
  std::vector<double> charge(1); plmd.cmd("setCharges",&charge[0]);
  std::vector<Vector> pos(1); plmd.cmd("setPositions",&pos[0][0]);
  std::vector<Vector> forces(1); plmd.cmd("setForces",&forces[0][0]);

  FILE* fp=fopen(ofile.c_str(),"r");
  bool do_read=true; unsigned nfram=0;
  std::vector<double> alig, disp; 
  while (do_read ) {
    PDB mypdb;
      // Read the pdb file
      do_read=mypdb.readFromFilepointer(fp,false,0.1);
      if( do_read ) {
          std::string num; Tools::convert( nfram+1, num ); nfram++;
          std::string iinput = "ref_" + num + ": READ_CONFIG REFERENCE=" + ofile + " NUMBER=" + num;
          if( argstr.length()>0 ) iinput += " READ_ARG=" + argstr;
          const char* icinp=iinput.c_str(); plmd.cmd("readInputLine",icinp);
          if( nfram==1 ) {
              alig.resize( mypdb.getOccupancy().size() );
              for(unsigned i=0;i<alig.size();++i) alig[i]=mypdb.getOccupancy()[i];
              disp.resize( mypdb.getBeta().size() );
              for(unsigned i=0;i<disp.size();++i) disp[i]=mypdb.getBeta()[i];
          }
      }
    }
    // Now create the objects to measure the distances between the frames
    std::vector<double> data( nfram );
    for(unsigned j=1;j<nfram;++j) {
        std::string istr, jstr; Tools::convert( j, istr); Tools::convert( j+1, jstr );
        if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) {
            std::string dstring = "d" + istr + ": RMSD TYPE=" + mtype + " REFERENCE_ATOMS=ref_" + istr +  " ATOMS=ref_" + jstr;
            // Get the align values 
            std::string anum; Tools::convert( alig[0], anum ); dstring += " ALIGN=" + anum;
            for(unsigned i=1;i<alig.size();++i){ Tools::convert( alig[i], anum ); dstring += "," + anum; }
            // Get the displace values
            std::string dnum; Tools::convert( disp[0], dnum ); dstring += " DISPLACE=" + dnum;
            for(unsigned i=1;i<disp.size();++i){ Tools::convert( disp[i], dnum ); dstring += "," + dnum; }
            const char* icinp=dstring.c_str(); plmd.cmd("readInputLine",icinp);
        } else if( mtype=="EUCLIDEAN" ) {
            // Difference between args
            std::string diff_string = "diff" + istr + ": DIFFERENCE ARG1=ref_" + istr +  " ARG2=ref_" + jstr;
            const char* icinp=diff_string.c_str(); plmd.cmd("readInputLine",icinp);
            // Retrieve number of quantities so we can make the powers vector
            ActionWithValue* aval = plmd.getActionSet().selectWithLabel<ActionWithValue*>( "diff" + istr );
            unsigned nquantities = aval->copyOutput(0)->getShape()[0];
            // SUM OF SQUARED
            std::string comb_string = "comb" + istr + ": COMBINE ARG=diff" + istr + " PERIODIC=NO POWERS=2";
            for(unsigned i=1;i<nquantities;++i) comb_string += ",2";
            const char* icinp2=comb_string.c_str(); plmd.cmd("readInputLine",icinp2);
            // SQUARE ROOT
            std::string fstr = "d" + istr + ": MATHEVAL ARG=comb" + istr + " FUNC=sqrt(x) PERIODIC=NO";
            const char* icinp3=fstr.c_str(); plmd.cmd("readInputLine",icinp3); 
        } else {
            plumed_merror("metric type " + mtype + " has not been implemented");
        }
        long rank; plmd.cmd("getDataRank d" + istr, &rank );
        if( rank!=0 ) error("distance should be of rank 1");
        std::vector<long> ishape(1); plmd.cmd("getDataShape d" + istr, &ishape[0] );  
        plmd.cmd("setMemoryForData d" + istr, &data[j-1] );
    }
    plmd.cmd("calc"); double mean=0;
    printf("N.B. THIS CODE ALWAYS AIMS TO CREATE EQUALLY SPACED FRAMES");
    printf("THERE MAY BE SMALL DESCREPENCIES IN THE NUMBERS BELOW, HOWEVER, BECAUSE OF ROUNDING ERRORS");
    for(unsigned i=1; i<nfram; ++i) {
        printf("FINAL DISTANCE BETWEEN FRAME %u AND %u IS %f \n",i-1,i,data[i-1]); mean += data[i-1];
    }
    printf("SUGGESTED LAMBDA PARAMETER IS THUS %f \n",2.3/mean/static_cast<double>( nfram-1 ) );
}

} // End of namespace
}
