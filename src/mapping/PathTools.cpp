/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "core/Value.h"
#include "reference/ReferenceConfiguration.h"
#include "PathReparameterization.h"
#include "reference/MetricRegister.h"
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
  std::string mtype; parse("--metric",mtype);
  std::string ifilename; parse("--path",ifilename);
  std::string ofmt; parse("--arg-fmt",ofmt);
  std::string ofilename; parse("--out",ofilename);
  if( ifilename.length()>0 ) {
    fprintf(out,"Reparameterising path in file named %s so that all frames are equally spaced \n",ifilename.c_str() );
    FILE* fp=fopen(ifilename.c_str(),"r");
    bool do_read=true; std::vector<std::unique_ptr<ReferenceConfiguration>> frames;
    while (do_read) {
      PDB mypdb;
      // Read the pdb file
      do_read=mypdb.readFromFilepointer(fp,false,0.1);
      if( do_read ) {
        auto mymsd(metricRegister().create<ReferenceConfiguration>( mtype, mypdb ));
        frames.emplace_back( std::move(mymsd) );
      }
    }
    std::vector<unsigned> fixed; parseVector("--fixed",fixed);
    if( fixed.size()==1 ) {
      if( fixed[0]!=0 ) error("input to --fixed should be two integers");
      fixed.resize(2); fixed[0]=0; fixed[1]=frames.size()-1;
    } else if( fixed.size()==2 ) {
      if( fixed[0]<0 || fixed[1]<0 || fixed[0]>(frames.size()-1) || fixed[1]>(frames.size()-1) ) {
        error("input to --fixed should be two numbers between 0 and the number of frames-1");
      }
    } else {
      error("input to --fixed should be two integers");
    }
    std::vector<AtomNumber> atoms; std::vector<std::string> arg_names;
    for(unsigned i=0; i<frames.size(); ++i) {
      frames[i]->getAtomRequests( atoms);
      frames[i]->getArgumentRequests( arg_names );
    }
    // Generate stuff to reparameterize
    Pbc fake_pbc; std::vector<std::unique_ptr<Value>> vals;
    for(unsigned i=0; i<frames[0]->getNumberOfReferenceArguments(); ++i) {
      vals.emplace_back(new Value()); vals[vals.size()-1]->setNotPeriodic();
    }

    // temporary pointes used to make the conversion once

    auto vals_ptr=Tools::unique2raw(vals);
    // And reparameterize
    PathReparameterization myparam( fake_pbc, vals_ptr, frames );
    // And make all points equally spaced
    double tol; parse("--tolerance",tol); myparam.reparameterize( fixed[0], fixed[1], tol );

    // Ouput data on spacings
    double mean=0;
    MultiValue myvpack( 1, frames[0]->getNumberOfReferenceArguments() + 3*frames[0]->getNumberOfReferencePositions() + 9 );
    ReferenceValuePack mypack( frames[0]->getNumberOfReferenceArguments(), frames[0]->getNumberOfReferencePositions(), myvpack );
    for(unsigned i=1; i<frames.size(); ++i) {
      double len = frames[i]->calc( frames[i-1]->getReferencePositions(), fake_pbc, vals_ptr, frames[i-1]->getReferenceArguments(), mypack, false );
      printf("FINAL DISTANCE BETWEEN FRAME %u AND %u IS %f \n",i-1,i,len );
      mean+=len;
    }
    printf("SUGGESTED LAMBDA PARAMETER IS THUS %f \n",2.3/mean/static_cast<double>( frames.size()-1 ) );

    // Delete all the frames
    OFile ofile; ofile.open(ofilename);
    std::vector<std::string> argnames; frames[0]->getArgumentRequests( argnames );
    std::vector<AtomNumber> atindices; frames[0]->getAtomRequests( atindices );
    PDB mypdb; mypdb.setAtomNumbers( atindices ); mypdb.setArgumentNames( argnames );
    for(unsigned i=0; i<frames.size(); ++i) {
      mypdb.setAtomPositions( frames[i]->getReferencePositions() );
      for(unsigned j=0; j<argnames.size(); ++j) mypdb.setArgumentValue( argnames[j], frames[i]->getReferenceArguments()[j] );
      ofile.printf("REMARK TYPE=%s\n",mtype.c_str() );
      mypdb.print( 10, NULL, ofile, ofmt );
    }
    // Delete the vals as we don't need them
    // for(unsigned i=0; i<vals.size(); ++i) delete vals[i];
    // Return as we are done
    return 0;
  }

// Read initial frame
  std::string istart; parse("--start",istart); FILE* fp2=fopen(istart.c_str(),"r"); PDB mystartpdb;
  if( istart.length()==0 ) error("input is missing use --istart + --iend or --path");
  if( !mystartpdb.readFromFilepointer(fp2,false,0.1) ) error("could not read fila " + istart);
  auto sframe=metricRegister().create<ReferenceConfiguration>( mtype, mystartpdb );
  fclose(fp2);

// Read final frame
  std::string iend; parse("--end",iend); FILE* fp1=fopen(iend.c_str(),"r"); PDB myendpdb;
  if( iend.length()==0 ) error("input is missing using --istart + --iend or --path");
  if( !myendpdb.readFromFilepointer(fp1,false,0.1) ) error("could not read fila " + iend);
  auto eframe=metricRegister().create<ReferenceConfiguration>( mtype, myendpdb );
  fclose(fp1);
// Get atoms and arg requests
  std::vector<AtomNumber> atoms; std::vector<std::string> arg_names;
  sframe->getAtomRequests( atoms); eframe->getAtomRequests( atoms);
  sframe->getArgumentRequests( arg_names ); eframe->getArgumentRequests( arg_names );

// Now read in the rest of the instructions
  unsigned nbefore, nbetween, nafter;
  parse("--nframes-before-start",nbefore); parse("--nframes",nbetween); parse("--nframes-after-end",nafter);
  nbetween++;
  fprintf(out,"Generating linear path connecting structure in file named %s to structure in file named %s \n",istart.c_str(),iend.c_str() );
  fprintf(out,"A path consisting of %u equally-spaced frames before the initial structure, %u frames between the intial and final structures "
          "and %u frames after the final structure will be created \n",nbefore,nbetween,nafter);

// Create a vector of arguments to use for calculating displacements
  Pbc fpbc;
  std::vector<std::unique_ptr<Value>> args;
  for(unsigned i=0; i<eframe->getNumberOfReferenceArguments(); ++i) {
    args.emplace_back(new Value()); args[args.size()-1]->setNotPeriodic();
  }

  // convert pointer once:
  auto args_ptr=Tools::unique2raw(args);

// Calculate the distance between the start and the end
  MultiValue myvpack( 1, sframe->getNumberOfReferenceArguments() + 3*sframe->getNumberOfReferencePositions() + 9);
  ReferenceValuePack mypack( sframe->getNumberOfReferenceArguments(), sframe->getNumberOfReferencePositions(), myvpack );
  double pathlen = sframe->calc( eframe->getReferencePositions(), fpbc, args_ptr, eframe->getReferenceArguments(), mypack, false );
// And the spacing between frames
  double delr = 1.0 / static_cast<double>( nbetween );
// Calculate the vector connecting the start to the end
  PDB mypdb; mypdb.setAtomNumbers( sframe->getAbsoluteIndexes() ); mypdb.addBlockEnd( sframe->getAbsoluteIndexes().size() );
  if( sframe->getArgumentNames().size()>0 ) mypdb.setArgumentNames( sframe->getArgumentNames() );
  Direction mydir(ReferenceConfigurationOptions("DIRECTION")); sframe->setupPCAStorage( mypack ); mydir.read( mypdb ); mydir.zeroDirection();
  sframe->extractDisplacementVector( eframe->getReferencePositions(), args_ptr, eframe->getReferenceArguments(), false, mydir );

// Now create frames
  OFile ofile; ofile.open(ofilename); unsigned nframes=0;
  Direction pos(ReferenceConfigurationOptions("DIRECTION")); pos.read( mypdb );
  for(int i=0; i<nbefore; ++i) {
    pos.setDirection( sframe->getReferencePositions(), sframe->getReferenceArguments() );
    pos.displaceReferenceConfiguration( -i*delr, mydir );
    mypdb.setAtomPositions( pos.getReferencePositions() );
    for(unsigned j=0; j<pos.getReferenceArguments().size(); ++j) mypdb.setArgumentValue( sframe->getArgumentNames()[j], pos.getReferenceArgument(j) );
    ofile.printf("REMARK TYPE=%s\n",mtype.c_str() );
    mypdb.print( 10, NULL, ofile, ofmt ); nframes++;
  }
  for(unsigned i=1; i<nbetween; ++i) {
    pos.setDirection( sframe->getReferencePositions(), sframe->getReferenceArguments() );
    pos.displaceReferenceConfiguration( i*delr, mydir );
    mypdb.setAtomPositions( pos.getReferencePositions() );
    for(unsigned j=0; j<pos.getReferenceArguments().size(); ++j) mypdb.setArgumentValue( sframe->getArgumentNames()[j], pos.getReferenceArgument(j) );
    ofile.printf("REMARK TYPE=%s\n",mtype.c_str() );
    mypdb.print( 10, NULL, ofile, ofmt ); nframes++;
  }
  for(unsigned i=0; i<nafter; ++i) {
    pos.setDirection( eframe->getReferencePositions(), eframe->getReferenceArguments() );
    pos.displaceReferenceConfiguration( i*delr, mydir );
    mypdb.setAtomPositions( pos.getReferencePositions() );
    for(unsigned j=0; j<pos.getReferenceArguments().size(); ++j) mypdb.setArgumentValue( sframe->getArgumentNames()[j], pos.getReferenceArgument(j) );
    ofile.printf("REMARK TYPE=%s\n",mtype.c_str() );
    mypdb.print( 10, NULL, ofile, ofmt ); nframes++;
  }

// double mean=0; printf("DISTANCE BETWEEN ORIGINAL FRAMES %f \n",pathlen);
// for(unsigned i=1;i<final_path.size();++i){
//    double len = final_path[i]->calc( final_path[i-1]->getReferencePositions(), fpbc, args, final_path[i-1]->getReferenceArguments(), mypack, false );
//    printf("FINAL DISTANCE BETWEEN FRAME %u AND %u IS %f \n",i-1,i,len );
//    mean+=len;
// }
// printf("SUGGESTED LAMBDA PARAMETER IS THUS %f \n",2.3/mean/static_cast<double>( final_path.size()-1 ) );

// Delete the args as we don't need them anymore
//  for(unsigned i=0; i<args.size(); ++i) delete args[i];
  ofile.close(); return 0;
}

} // End of namespace
}
