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
#include "Mapping.h"
#include "TrigonometricPathVessel.h"
#include "PathReparameterization.h"
#include "reference/Direction.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/SetupMolInfo.h"

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

The input below provides an example that shows how the adaptive path works. The path is updated every 50 steps of
MD based on the data accumulated during the preceding 50 time steps.

\plumedfile
d1: DISTANCE ATOMS=1,2 COMPONENTS
pp: ADAPTIVE_PATH TYPE=EUCLIDEAN FIXED=2,5 UPDATE=50 WFILE=out-path.pdb WSTRIDE=50 REFERENCE=mypath.pdb
PRINT ARG=d1.x,d1.y,pp.* FILE=colvar
\endplumedfile

In the case above the distance between frames is calculated based on the \f$x\f$ and \f$y\f$ components of the vector connecting
atoms 1 and 2.  As such an extract from the input reference path (mypath.pdb) would look as follows:

\auxfile{mypath.pdb}
REMARK ARG=d1.x,d1.y d1.x=1.12 d1.y=-.60
END
REMARK ARG=d1.x,d1.y d1.x=.99 d1.y=-.45
END
REMARK ARG=d1.x,d1.y d1.x=.86 d1.y=-.30
END
REMARK ARG=d1.x,d1.y d1.x=.73 d1.y=-.15
END
REMARK ARG=d1.x,d1.y d1.x=.60 d1.y=0
END
REMARK ARG=d1.x,d1.y d1.x=.47 d1.y=.15
END
\endauxfile

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
  TrigonometricPathVessel* mypathv;
  std::vector<double> wsum;
  Direction displacement,displacement2;
  std::vector<Direction> pdisplacements;
public:
  static void registerKeywords( Keywords& keys );
  explicit AdaptivePath(const ActionOptions&);
  void calculate() override;
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const override;
  double getLambda() override { return 0.0; }
  double transformHD( const double& dist, double& df ) const override;
  void update() override;
};

PLUMED_REGISTER_ACTION(AdaptivePath,"ADAPTIVE_PATH")

void AdaptivePath::registerKeywords( Keywords& keys ) {
  Mapping::registerKeywords( keys ); keys.remove("PROPERTY");
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
  setLowMemOption( true ); parseVector("FIXED",fixedn);
  if( fixedn[0]<1 || fixedn[1]>getNumberOfReferencePoints() ) error("fixed nodes must be in range from 0 to number of nodes");
  if( fixedn[0]>=fixedn[1] ) error("invalid selection for fixed nodes first index provided must be smaller than second index");
  log.printf("  fixing position of frames numbered %u and %u \n",fixedn[0],fixedn[1]);
  fixedn[0]--; fixedn[1]--;   // Set fixed notes with c++ indexing starting from zero
  parse("UPDATE",update_str); if( update_str<1 ) error("update frequency for path should be greater than or equal to one");
  log.printf("  updating path every %u MD steps \n",update_str);

  double halflife; parse("HALFLIFE",halflife);
  if( halflife<0 ) fadefact=1.0;
  else {
    fadefact = exp( -0.693147180559945 / static_cast<double>(halflife) );
    log.printf("  weight of contribution to frame halves every %f steps \n",halflife);
  }

  // Create the list of tasks (and reset projections of frames)
  PDB mypdb; mypdb.setAtomNumbers( getAbsoluteIndexes() ); mypdb.addBlockEnd( getAbsoluteIndexes().size() );
  std::vector<std::string> argument_names( getNumberOfArguments() );
  for(unsigned i=0; i<getNumberOfArguments(); ++i) argument_names[i] = getPntrToArgument(i)->getName();
  if( argument_names.size()>0 ) mypdb.setArgumentNames( argument_names );
  displacement.read( mypdb ); displacement2.read( mypdb );
  for(int i=0; i<getNumberOfReferencePoints(); ++i) {
    addTaskToList( i ); pdisplacements.push_back( Direction(ReferenceConfigurationOptions("DIRECTION")) );
    property.find("spath")->second[i] = static_cast<double>( i - static_cast<int>(fixedn[0]) ) / static_cast<double>( fixedn[1] - fixedn[0] );
    pdisplacements[i].read( mypdb ); wsum.push_back( 0.0 );
  }
  plumed_assert( getPropertyValue( fixedn[0], "spath" )==0.0 && getPropertyValue( fixedn[1], "spath" )==1.0 );
  // And activate them all
  deactivateAllTasks();
  for(unsigned i=0; i<getFullNumberOfTasks(); ++i) taskFlags[i]=1;
  lockContributors();

  // Setup the vessel to hold the trig path
  std::string input; addVessel("GPATH", input, -1 );
  readVesselKeywords();
  // Check that there is only one vessel - the one holding the trig path
  plumed_dbg_assert( getNumberOfVessels()==1 );
  // Retrieve the path vessel
  mypathv = dynamic_cast<TrigonometricPathVessel*>( getPntrToVessel(0) );
  plumed_assert( mypathv );

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

void AdaptivePath::calculate() {
  runAllTasks();
}

void AdaptivePath::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  // This builds a pack to hold the derivatives
  ReferenceValuePack mypack( getNumberOfArguments(), getNumberOfAtoms(), myvals );
  finishPackSetup( current, mypack );
  // Calculate the distance from the frame
  double val=calculateDistanceFunction( current, mypack, true );
  // Put the element value in element zero
  myvals.setValue( 0, val ); myvals.setValue( 1, 1.0 );
  return;
}

double AdaptivePath::transformHD( const double& dist, double& df ) const {
  df=1.0; return dist;
}

void AdaptivePath::update() {
  double weight2 = -1.*mypathv->dx;
  double weight1 = 1.0 + mypathv->dx;
  if( weight1>1.0 ) {
    weight1=1.0; weight2=0.0;
  } else if( weight2>1.0 ) {
    weight1=0.0; weight2=1.0;
  }
  // Add projections to dispalcement accumulators
  ReferenceConfiguration* myref = getReferenceConfiguration( mypathv->iclose1 );
  myref->extractDisplacementVector( getPositions(), getArguments(), mypathv->cargs, false, displacement );
  getReferenceConfiguration( mypathv->iclose2 )->extractDisplacementVector( myref->getReferencePositions(), getArguments(), myref->getReferenceArguments(), false, displacement2 );
  displacement.addDirection( -mypathv->dx, displacement2 );
  pdisplacements[mypathv->iclose1].addDirection( weight1, displacement );
  pdisplacements[mypathv->iclose2].addDirection( weight2, displacement );
  // Update weight accumulators
  wsum[mypathv->iclose1] *= fadefact;
  wsum[mypathv->iclose2] *= fadefact;
  wsum[mypathv->iclose1] += weight1;
  wsum[mypathv->iclose2] += weight2;

  // This does the update of the path if it is time to
  if( (getStep()>0) && (getStep()%update_str==0) ) {
    wsum[fixedn[0]]=wsum[fixedn[1]]=0.;
    for(unsigned inode=0; inode<getNumberOfReferencePoints(); ++inode) {
      if( wsum[inode]>0 ) {
        // First displace the node by the weighted direction
        getReferenceConfiguration( inode )->displaceReferenceConfiguration( 1./wsum[inode], pdisplacements[inode] );
        // Reset the displacement
        pdisplacements[inode].zeroDirection();
      }
    }
    // Now ensure all the nodes of the path are equally spaced
    PathReparameterization myspacings( getPbc(), getArguments(), getAllReferenceConfigurations() );
    myspacings.reparameterize( fixedn[0], fixedn[1], tolerance );
  }
  if( (getStep()>0) && (getStep()%wstride==0) ) {
    pathfile<<"# PATH AT STEP "<<getStep();
    pathfile.printf(" TIME %f \n",getTime());
    std::vector<std::unique_ptr<ReferenceConfiguration>>& myconfs=getAllReferenceConfigurations();
    std::vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
    if( moldat.size()>1 ) error("you should only have one MOLINFO action in your input file");
    SetupMolInfo* mymoldat=NULL; if( moldat.size()==1 ) mymoldat=moldat[0];
    std::vector<std::string> argument_names( getNumberOfArguments() );
    for(unsigned i=0; i<getNumberOfArguments(); ++i) argument_names[i] = getPntrToArgument(i)->getName();
    PDB mypdb; mypdb.setArgumentNames( argument_names );
    for(unsigned i=0; i<myconfs.size(); ++i) {
      pathfile.printf("REMARK TYPE=%s\n", myconfs[i]->getName().c_str() );
      mypdb.setAtomPositions( myconfs[i]->getReferencePositions() );
      for(unsigned j=0; j<getNumberOfArguments(); ++j) mypdb.setArgumentValue( getPntrToArgument(j)->getName(), myconfs[i]->getReferenceArgument(j) );
      mypdb.print( atoms.getUnits().getLength()/0.1, mymoldat, pathfile, ofmt );
    }
    pathfile.flush();
  }
}

}
}
