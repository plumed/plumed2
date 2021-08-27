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
#include "core/ActionShortcut.h"
#include "core/AverageBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace dimred {

//+PLUMEDOC DIMRED PCA
/*
Perform principal component analysis (PCA) using either the positions of the atoms a large number of collective variables as input.

Principal component analysis is a statistical technique that uses an orthogonal transformation to convert a set of observations of
poorly correlated variables into a set of linearly uncorrelated variables.  You can read more about the specifics of this technique
here: https://en.wikipedia.org/wiki/Principal_component_analysis

When used with molecular dynamics simulations a set of frames taken from the trajectory, \f$\{X_i\}\f$, or the values of
a number of collective variables which are calculated from the trajectory frames are used as input.  In this second instance your
input to the PCA analysis algorithm is thus a set of high-dimensional vectors of collective variables.  However, if
collective variables are calculated from the positions of the atoms or if the positions are used directly the assumption is that
this input trajectory is a set of poorly correlated (high-dimensional) vectors.  After principal component analysis has been
performed the output is a set of orthogonal vectors that describe the directions in which the largest motions have been seen.
In other words, principal component analysis provides a method for lowering the dimensionality of the data contained in a trajectory.
These output directions are some linear combination of the \f$x\f$, \f$y\f$ and \f$z\f$ positions if the positions were used as input
or some linear combination of the input collective variables if a high-dimensional vector of collective variables was used as input.

As explained on the Wikipedia page you must calculate the average and covariance for each of the input coordinates.  In other words, you must
calculate the average structure and the amount the system fluctuates around this average structure.  The problem in doing so when the
\f$x\f$, \f$y\f$ and \f$z\f$ coordinates of a molecule are used as input is that the majority of the changes in the positions of the
atoms comes from the translational and rotational degrees of freedom of the molecule.  The first six principal components will thus, most likely,
be uninteresting.  Consequently, to remedy this problem PLUMED provides the functionality to perform an RMSD alignment of the all the structures
to be analyzed to the first frame in the trajectory.  This can be used to effectively remove translational and/or rotational motions from
consideration.  The resulting principal components thus describe vibrational motions of the molecule.

If you wish to calculate the projection of a trajectory on a set of principal components calculated from this PCA action then the output can be
used as input for the \ref PCAVARS action.

\par Examples

The following input instructs PLUMED to perform a principal component analysis in which the covariance matrix is calculated from changes in the positions
of the first 22 atoms.  The TYPE=OPTIMAL instruction ensures that translational and rotational degrees of freedom are removed from consideration.
The first two principal components will be output to a file called PCA-comp.pdb.  Trajectory frames will be collected on every step and the PCA calculation
will be performed at the end of the simulation.

\plumedfile
ff: COLLECT_FRAMES ATOMS=1-22 STRIDE=1
pca: PCA USE_OUTPUT_DATA_FROM=ff METRIC=OPTIMAL NLOW_DIM=2
OUTPUT_PCA_PROJECTION USE_OUTPUT_DATA_FROM=pca FILE=PCA-comp.pdb
\endplumedfile

The following input instructs PLUMED to perform a principal component analysis in which the covariance matrix is calculated from changes in the six distances
seen in the previous lines.  Notice that here the TYPE=EUCLIDEAN keyword is used to indicate that no alignment has to be done when calculating the various
elements of the covariance matrix from the input vectors.  In this calculation the first two principal components will be output to a file called PCA-comp.pdb.
Trajectory frames will be collected every five steps and the PCA calculation is performed every 1000 steps.  Consequently, if you run a 2000 step simulation the
PCA analysis will be performed twice.  The REWEIGHT_BIAS action in this input tells PLUMED that rather that ascribing a weight of one to each of the frames
when calculating averages and covariance matrices a reweighting should be performed based and each frames' weight in these calculations should be determined based on
the current value of the instantaneous bias (see \ref REWEIGHT_BIAS).

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,3
d3: DISTANCE ATOMS=1,4
d4: DISTANCE ATOMS=2,3
d5: DISTANCE ATOMS=2,4
d6: DISTANCE ATOMS=3,4
rr: RESTRAINT ARG=d1 AT=0.1 KAPPA=10
rbias: REWEIGHT_BIAS TEMP=300

ff: COLLECT_FRAMES ARG=d1,d2,d3,d4,d5,d6 LOGWEIGHTS=rbias STRIDE=5
pca: PCA USE_OUTPUT_DATA_FROM=ff METRIC=EUCLIDEAN NLOW_DIM=2
OUTPUT_PCA_PROJECTION USE_OUTPUT_DATA_FROM=pca STRIDE=100 FILE=PCA-comp.pdb
\endplumedfile

*/
//+ENDPLUMEDOC

class PCA : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  PCA( const ActionOptions& );
};

PLUMED_REGISTER_ACTION(PCA,"PCA")

void PCA::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the arguments that you would like to make the histogram for");
  keys.add("compulsory","NLOW_DIM","number of low-dimensional coordinates required"); 
  keys.add("optional","FILE","the file on which to output the low dimensional coordinates");
  keys.add("optional","FMT","the format to use when outputting the low dimensional coordinates");
}


PCA::PCA( const ActionOptions& ao ) :
  Action(ao),
  ActionShortcut(ao)
{
  // Read in the data set we are doing PCA for
  std::string arg; parse("ARG",arg);
  AverageBase* mydata = plumed.getActionSet().selectWithLabel<AverageBase*>(arg);
  if( !mydata ) error("input to PCA should be a COLLECT_FRAMES or COLLECT_REPLICAS object");
  std::string argstr; unsigned anum=1; bool hasargs=false;
  for(unsigned i=0;i<mydata->getNumberOfComponents();++i) {
      std::string thislab = mydata->copyOutput(i)->getName();
      if( thislab.find(".logweights")==std::string::npos && thislab.find(".pos")==std::string::npos ) {
          std::string num; Tools::convert( anum, num ); 
          // Average for component
          std::size_t dot = thislab.find("."); std::string datalab = thislab.substr(dot+1);
          // This calculates the average
          readInputLine( getShortcutLabel() + "_average" + num + ": AVERAGE ARG=" + datalab + mydata->getStrideClearAndWeights() );
          // This calculates the centered data vector
          readInputLine( getShortcutLabel() + "_centeredvec" + num + ": MATHEVAL ARG1=" + thislab + " ARG2=" + getShortcutLabel() + "_average" + num + " FUNC=x-y PERIODIC=NO");
          argstr += " ARG" + num + "=" + getShortcutLabel() + "_centeredvec" + num; hasargs=true; anum++; 
      }
  }
  // Additional lines will be added here to deal with the atoms
  if( mydata->getNumberOfAtoms()>0 ) {
      // This calculates the average position
      readInputLine( getShortcutLabel() + "_average_atoms: AVERAGE " + mydata->getAtomsData() + " " +  mydata->getStrideClearAndWeights() );
      for(unsigned i=0;i<mydata->getNumberOfAtoms();++i) {
          std::string num, pnum, atnum; Tools::convert( i+1, atnum );
          // This calculates the centered data vector for the x component
          Tools::convert( anum, num ); Tools::convert( 3*i + 1, pnum );
          readInputLine( getShortcutLabel() + "_avpos" + num + ": SELECT_COMPONENTS ARG=" + getShortcutLabel() + "_average_atoms COMPONENTS=" + pnum );
          readInputLine( getShortcutLabel() + "_centeredvec" + num + ": MATHEVAL ARG1=" + arg + ".posx-" + atnum + 
                                                                               " ARG2=" + getShortcutLabel() + "_avpos"  + pnum + 
                                                                               " FUNC=x-y PERIODIC=NO"); 
          argstr += " ARG" + num + "=" + getShortcutLabel() + "_centeredvec" + num; anum++;
          // This calculates the centered data vector for the y component 
          Tools::convert( anum, num ); Tools::convert( 3*i + 2, pnum );
          readInputLine( getShortcutLabel() + "_avpos" + num + ": SELECT_COMPONENTS ARG=" + getShortcutLabel() + "_average_atoms COMPONENTS=" + pnum ); 
          readInputLine( getShortcutLabel() + "_centeredvec" + num + ": MATHEVAL ARG1=" + arg + ".posy-" + atnum +
                                                                               " ARG2=" + getShortcutLabel() + "_avpos"  + pnum + 
                                                                               " FUNC=x-y PERIODIC=NO");
          argstr += " ARG" + num + "=" + getShortcutLabel() + "_centeredvec" + num; anum++;
          // This calculates the centered data vector for the z component 
          Tools::convert( anum, num ); Tools::convert( 3*i + 3, pnum ); 
          readInputLine( getShortcutLabel() + "_avpos" + num + ": SELECT_COMPONENTS ARG=" + getShortcutLabel() + "_average_atoms COMPONENTS=" + pnum );
          readInputLine( getShortcutLabel() + "_centeredvec" + num + ": MATHEVAL ARG1=" + arg + ".posz-" + atnum +
                                                                               " ARG2=" + getShortcutLabel() + "_avpos"  + pnum +
                                                                               " FUNC=x-y PERIODIC=NO");
          argstr += " ARG" + num + "=" + getShortcutLabel() + "_centeredvec" + num; anum++;
      } 
  }
  // Now antilog the weights
  readInputLine( getShortcutLabel() + "_weights: MATHEVAL ARG1=" + arg + ".logweights FUNC=exp(x) PERIODIC=NO");
  // And calculate the covariance matrix
  readInputLine( getShortcutLabel() + "_covar: COVARIANCE_MATRIX " + argstr + " WEIGHTS=" + getShortcutLabel() + "_weights");
  // Diagonalize the covariance matrix
  unsigned ndim; parse("NLOW_DIM",ndim); std::string vecstr="1";
  if( ndim<=0 || ndim>mydata->getNumberOfComponents()-1 ) error("cannot generate projection in space of dimension higher than input coordinates"); 
  for(unsigned i=1;i<ndim;++i){ std::string num; Tools::convert( i+1, num ); vecstr += "," + num; }
  readInputLine( getShortcutLabel() + "_eig: DIAGONALIZE ARG=" + getShortcutLabel() + "_covar VECTORS=" + vecstr );
  // Now set up printing 
  std::string filename; parse("FILE",filename);
  if( filename.length()>0 ) {
      if( filename.find(".pdb")==std::string::npos ) error("output file for PCA should be a PDB");
      std::string fmt; parse("FMT",fmt); if( fmt.length()==0 ) fmt="%f";
      std::string atstr, avlist;
      if( hasargs ) {
         avlist = " CONFIG1=" + getShortcutLabel() + "_average1"; 
         for(unsigned i=1;i<anum-1;++i) { std::string num; Tools::convert( i+1, num ); avlist += "," + getShortcutLabel() + "_average" + num; }
         if( mydata->getNumberOfAtoms()>0 ) avlist = "," + getShortcutLabel() + "_average_atoms";
      }
      if( mydata->getNumberOfAtoms()>0 ) atstr = " CONFIG1=" + getShortcutLabel() + "_average_atoms";
      std::string eiglist; 
      for(unsigned i=0;i<ndim;++i) {
          std::string lnum, num; Tools::convert( i+2, lnum ); Tools::convert( i+1, num ); 
          eiglist += " CONFIG" + lnum + "=" + getShortcutLabel() + "_eig.vecs-" + num;
      } 
      readInputLine("PRINT DESCRIPTION=PCA FILE=" + filename + " FMT=" + fmt + atstr + avlist + " " + eiglist );
  }
  // And calculate the projections of the stored data on to the PCA vectors
  for(unsigned i=0;i<ndim;++i) {
      std::string num; Tools::convert( i+1, num );
      readInputLine( getShortcutLabel() + "-" + num + ": PROJECT_ON_VECTOR " + argstr + " VECTOR=" + getShortcutLabel() + "_eig.vecs-" + num ); 
  }
}

}
}
