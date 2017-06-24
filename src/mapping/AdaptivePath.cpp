/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
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
#include "Mapping.h"
#include "PathReparameterization.h"
#include "reference/Direction.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

//+PLUMEDOC COLVAR ADAPTIVE_PATH
/*
Compute path collective variables that adapt to the lowest free energy path connecting states A and B.

The Path Collective Variables developed by Branduardi and co-workers \cite brand07 allow one
to compute the progress along a high-dimensional path and the distance from the high-dimensional
path.  The progress along the path (s) is computed using:

\f[
s = i_2 + \textrm{sign}(i_2-i_1) \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2}
\f]

In this expression \f$\mathbf{v}_1\f$ and \f$\mathbf{v}_3\f$ are the vectors connecting the current position to the closest and second closest node of the path,
respectfully and \f$i_1\f$ and \f$i_2\f$ are the projections of the closest and second closest frames of the path. \f$\mathbf{v}_2\f$, meanwhile, is the
vector connecting the closest frame to the second closest frame.  The distance from the path, \f$z\f$ is calculated using:

\f[
z = \sqrt{ \left[ |\mathbf{v}_1|^2 - |\mathbf{v}_2| \left( \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2} \right) \right]^2 }
\f]

Notice that these are the definitions of \f$s\f$ and \f$z\f$ that are used by \ref PATH when the GPATH option is employed.  The reason for this is that
the adaptive path method implemented in this action was inspired by the work of Diaz and Ensing in which these formula were used \cite BerndAdaptivePath.
To learn more about how the path is adapted we strongly recommend reading this paper.

\par Examples

The input below provides an example of how the adaptive path works in practise. The path is updated every 50 steps of
MD based on the data accumulated during the preceding 50 time steps.

\plumedfile
d1: DISTANCE ATOMS=1,2 COMPONENTS
pp: ADAPTIVE_PATH TYPE=EUCLIDEAN FIXED=5,15 UPDATE=50 WFILE=out-path.pdb WSTRIDE=50 REFERENCE=mypath.pdb
PRINT ARG=d1.x,d1.y,pp.* FILE=colvar
\endplumedfile

In the case above the distance between frames is calculated based on the \f$x\f$ and \f$y\f$ components of the vector connecting
atoms 1 and 2.  As such an extract from the input reference path (mypath.pdb) would look as follows:

\verbatim
REMARK ARG=d1.x,d1.y d1.x=1.12 d1.y=-.60
END
REMARK ARG=d1.x,d1.y d1.x=.99 d1.y=-.45
END
\endverbatim

Notice that one can also use RMSD frames in place of arguments like those above.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class AdaptivePath : public Mapping {
private:
  OFile pathfile;
  std::string ofmt;
  double fadefact, tolerance;
  unsigned update_str, wstride;
  std::vector<unsigned> fixedn;
  std::vector<double> wsum;
  Direction displacement, displacement2;
  std::vector<Direction> pdisplacements;
public:
  static void registerKeywords( Keywords& keys );
  explicit AdaptivePath(const ActionOptions&);
  void update();
};

PLUMED_REGISTER_ACTION(AdaptivePath,"ADAPTIVE_PATH")

void AdaptivePath::registerKeywords( Keywords& keys ) {
  Mapping::registerKeywords( keys ); keys.remove("SQUARED"); 
  keys.add("compulsory","FIXED","the positions in the list of input frames of the two path nodes whose positions remain fixed during the path optimization");
  keys.add("compulsory","HALFLIFE","-1","the number of MD steps after which a previously measured path distance weighs only 50% in the average. This option may increase convergence by allowing to \"forget\" the memory of a bad initial guess path. The default is to set this to infinity");
  keys.add("compulsory","UPDATE","the frequency with which the path should be updated");
  keys.add("compulsory","TOLERANCE","1E-6","the tolerance to use for the path updating algorithm that makes all frames equidistant");
  keys.add("optional","WFILE","file on which to write out the path");
  keys.add("compulsory","FMT","%f","the format to use for output files");
  keys.add("optional","WSTRIDE","frequency with which to write out the path");
}

AdaptivePath::AdaptivePath(const ActionOptions& ao):
  Action(ao),
  Mapping(ao),
  fixedn(2),
  displacement(ReferenceConfigurationOptions("DIRECTION")),
  displacement2(ReferenceConfigurationOptions("DIRECTION"))
{
  if( myframes.size()<2 ) error("not enough reference configurations in path");
  parseVector("FIXED",fixedn);
  if( fixedn[0]<1 || fixedn[1]>myframes.size() ) error("fixed nodes must be in range from 0 to number of nodes");
  if( fixedn[0]>=fixedn[1] ) error("invalid selection for fixed nodes first index provided must be smaller than second index");
  log.printf("  fixing position of frames numbered %u and %u \n",fixedn[0],fixedn[1]);
  fixedn[0]--; fixedn[1]--;   // Set fixed notes with c++ indexing starting from zero
  parse("UPDATE",update_str); if( update_str<1 ) error("update frequency for path should be greater than or equal to one");
  log.printf("  updating path every %u MD steps \n",update_str);

  double halflife; parse("HALFLIFE",halflife);
  log.printf("  weight of contribution to frame halves every %f steps \n",halflife);
  if( halflife<0 ) fadefact=1.0;
  else fadefact = exp( -0.693147180559945 / static_cast<double>(halflife) );

  // Create the list of tasks (and reset projections of frames)
  std::vector<std::string> argument_names( getNumberOfArguments() );
  for(unsigned i=0; i<getNumberOfArguments(); ++i) argument_names[i] = getPntrToArgument(i)->getName();
  displacement.setNamesAndAtomNumbers( getAbsoluteIndexes(), argument_names );
  displacement2.setNamesAndAtomNumbers( getAbsoluteIndexes(), argument_names );
  for(int i=0; i<myframes.size(); ++i) {
    pdisplacements.push_back( Direction(ReferenceConfigurationOptions("DIRECTION")) );
    pdisplacements[i].setNamesAndAtomNumbers( getAbsoluteIndexes(), argument_names ); wsum.push_back( 0.0 );
  }
  // Make sure we collect all the data
  squared = true; getPntrToComponent(0)->buildDataStore();

  // Information for write out
  std::string wfilename; parse("WFILE",wfilename);
  if( wfilename.length()>0 ) {
    wstride=0; parse("WSTRIDE",wstride); parse("FMT",ofmt);
    pathfile.link(*this); pathfile.open( wfilename ); pathfile.setHeavyFlush();
    if( wstride<update_str ) error("makes no sense to write out path more frequently than update stride");
    log.printf("  writing path out every %u steps to file named %s with format %s \n",wstride,wfilename.c_str(),ofmt.c_str());
  }
  log<<"  Bibliography "<<plumed.cite("Diaz Leines and Ensing, Phys. Rev. Lett. 109, 020601 (2012)")<<"\n";
}

void AdaptivePath::update() {
  runAllTasks();   // Redo calculation - makes sure no problems due to parallel tempering

  double v1v1 = getPntrToComponent(0)->get(0); unsigned iclose1 = 0;
  double v3v3 = getPntrToComponent(0)->get(1); unsigned iclose2 = 1;
  if( v1v1>v3v3 ){
      double tmp=v1v1; v1v1=v3v3; v3v3=tmp;
      iclose1 = 1; iclose2 = 0;
  }
  for(unsigned i=2; i<myframes.size(); ++i) {
    double ndist=getPntrToComponent(0)->get(i);
    if( ndist<v1v1 ) {
      v3v3=v1v1; iclose2=iclose1;
      v1v1=ndist; iclose1=i;
    } else if( ndist<v3v3 ) {
      v3v3=ndist; iclose2=i;
    }
  }
  // And find third closest point
  int isign = iclose1 - iclose2;
  if( isign>1 ) isign=1; else if( isign<-1 ) isign=-1;
  int iclose3 = iclose1 + isign; double v2v2;

  //  Create some holders for stuff
  MultiValue mydpack1(1,0), mydpack2(1,0); 
  ReferenceValuePack mypack1( getNumberOfArguments(), getNumberOfAtoms(), mydpack1 ); myframes[0]->setupPCAStorage( mypack1 );
  ReferenceValuePack mypack2( getNumberOfArguments(), getNumberOfAtoms(), mydpack2 );
  mypack1.clear(); calculateDistanceFromReference( iclose1, mypack1 );
  if( iclose3<0 || iclose3>=myframes.size() ) {
    v2v2=calculateDistanceBetweenReferenceAndThisPoint( iclose2, myframes[iclose1]->getReferencePositions(), myframes[iclose1]->getReferenceArguments(), mypack2 );
    extractDisplacementVector( iclose2, myframes[iclose1]->getReferencePositions(), myframes[iclose1]->getReferenceArguments(), displacement );
  } else {
    v2v2=calculateDistanceBetweenReferenceAndThisPoint( iclose1, myframes[iclose3]->getReferencePositions(), myframes[iclose3]->getReferenceArguments(), mypack2 );
    extractDisplacementVector( iclose1, myframes[iclose3]->getReferencePositions(), myframes[iclose3]->getReferenceArguments(), displacement );
  }
  if( getNumberOfAtoms()>0 ) {
      ReferenceAtoms* at = dynamic_cast<ReferenceAtoms*>( myframes[iclose1] );
      const std::vector<double> & displace( at->getDisplace() );
      for(unsigned i=0; i<getNumberOfAtoms(); ++i) mypack1.getAtomsDisplacementVector()[i] /= displace[i];
  }
  // Calculate the dot product of v1 with v2
  double v1v2 = projectDisplacementOnVector( iclose1, displacement, mypack1 );
  double root = sqrt( v1v2*v1v2 - v2v2 * ( v1v1 - v3v3) );
  double dx = 0.5 * ( (root + v1v2) / v2v2 - 1.);
  double weight2 = -1.* dx; double weight1 = 1.0 + dx;
  if( weight1>1.0 ) {
    weight1=1.0; weight2=0.0;
  } else if( weight2>1.0 ) {
    weight1=0.0; weight2=1.0;
  }
  // Add projections to dispalcement accumulators
  //ReferenceConfiguration* myref = getReferenceConfiguration( mypathv->iclose1 );
  // myframes[iclose1]->
  std::vector<double> cargs( getNumberOfScalarArguments() ); 
  for(unsigned i=0;i<getNumberOfScalarArguments();++i) cargs[i]=getArgumentScalar(i);
  extractDisplacementVector( iclose1, getPositions(), cargs, displacement );
  extractDisplacementVector( iclose2, myframes[iclose1]->getReferencePositions(), myframes[iclose1]->getReferenceArguments(), displacement2 );
//   myref->extractDisplacementVector( getPositions(), getArguments(), mypathv->cargs, false, displacement );
//   getReferenceConfiguration( iclose2 )->extractDisplacementVector( myref->getReferencePositions(), getArguments(), myref->getReferenceArguments(), false, displacement2 );
  displacement.addDirection( -dx, displacement2 );
  pdisplacements[iclose1].addDirection( weight1, displacement );
  pdisplacements[iclose2].addDirection( weight2, displacement );
  // Update weight accumulators
  wsum[iclose1] *= fadefact;
  wsum[iclose2] *= fadefact;
  wsum[iclose1] += weight1;
  wsum[iclose2] += weight2;

  // This does the update of the path if it is time to
  if( (getStep()>0) && (getStep()%update_str==0) ) {
    wsum[fixedn[0]]=wsum[fixedn[1]]=0.;
    for(unsigned inode=0; inode<myframes.size(); ++inode) {
      if( wsum[inode]>0 ) {
        // First displace the node by the weighted direction
        myframes[inode]->displaceReferenceConfiguration( 1./wsum[inode], pdisplacements[inode] );
        // Reset the displacement
        pdisplacements[inode].zeroDirection();
      }
    }
    // Now ensure all the nodes of the path are equally spaced
    PathReparameterization myspacings( getPbc(), getArguments(), myframes );
    myspacings.reparameterize( fixedn[0], fixedn[1], tolerance );
  }
  if( (getStep()>0) && (getStep()%wstride==0) ) {
    pathfile.printf("# PATH AT STEP %d TIME %f \n", getStep(), getTime() );
    for(unsigned i=0; i<myframes.size(); ++i) myframes[i]->print( pathfile, ofmt, atoms.getUnits().getLength()/0.1 );
    pathfile.flush();
  }
}

}
}
