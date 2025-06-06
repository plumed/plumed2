/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "core/ActionAtomistic.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC DIMRED PCA
/*
Perform principal component analysis (PCA) using either the positions of the atoms a large number of collective variables as input.

Principal component analysis is a statistical technique that uses an orthogonal transformation to convert a set of observations of
poorly correlated variables into a set of linearly uncorrelated variables.  You can read more about the specifics of this technique
[here](https://en.wikipedia.org/wiki/Principal_component_analysis)

When used with molecular dynamics simulations a set of frames taken from the trajectory, $\{X_i\}$, or the values of
a number of collective variables which are calculated from the trajectory frames are used as input.  In this second instance your
input to the PCA analysis algorithm is thus a set of high-dimensional vectors of collective variables.  However, if
collective variables are calculated from the positions of the atoms or if the positions are used directly the assumption is that
this input trajectory is a set of poorly correlated (high-dimensional) vectors.  After principal component analysis has been
performed the output is a set of orthogonal vectors that describe the directions in which the largest motions have been seen.
In other words, principal component analysis provides a method for lowering the dimensionality of the data contained in a trajectory.
These output directions are some linear combination of the $x$, $y$ and $z$ positions if the positions were used as input
or some linear combination of the input collective variables if a high-dimensional vector of collective variables was used as input.

As explained on the Wikipedia page you must calculate the average and covariance for each of the input coordinates.  In other words, you must
calculate the average structure and the amount the system fluctuates around this average structure.  The problem in doing so when the
$x$, $y$ and $z$ coordinates of a molecule are used as input is that the majority of the changes in the positions of the
atoms comes from the translational and rotational degrees of freedom of the molecule.  The first six principal components will thus, most likely,
be uninteresting.  Consequently, to remedy this problem PLUMED provides the functionality to perform an RMSD alignment of the all the structures
to be analyzed to the first frame in the trajectory.  This can be used to effectively remove translational and/or rotational motions from
consideration.  The resulting principal components thus describe vibrational motions of the molecule.

If you wish to calculate the projection of a trajectory on a set of principal components calculated from this PCA action then the output can be
used as input for the [PCAVARS](PCAVARS.md) action.

## Examples

The following input instructs PLUMED to perform a principal component analysis in which the covariance matrix is calculated from changes in the positions
of the first 22 atoms. The `TYPE=OPTIMAL` instruction ensures that translational and rotational degrees of freedom are removed from consideration.
The average position and the first two principal components will be output to a file called `pca-comp.pdb`. Trajectory frames will be collected on every step and the PCA calculation
will be performed at the end of the simulation. The `colvar` file that is output contains the projections of all the positions in the high dimensional space on these vectors.

```plumed
ff: COLLECT_FRAMES ATOMS=1-22 STRIDE=1
pca: PCA ARG=ff NLOW_DIM=2 FILE=pca-comp.pdb
DUMPVECTOR ARG=pca,pca_weights FILE=colvar STRIDE=0
```

The following input instructs PLUMED to perform a principal component analysis in which the covariance matrix is calculated from changes in the six distances
seen in the input file below.  In this calculation the first two principal components will be output to a file called PCA-comp.pdb.
Trajectory frames will be collected every five steps and the PCA calculation is performed every 1000 steps. Consequently, if you run a 2000 step simulation the
PCA analysis will be performed twice. The REWEIGHT_BIAS action in this input tells PLUMED that rather that ascribing a weight of one to each of the frames
when calculating averages and covariance matrices a reweighting should be performed based and each frames' weight in these calculations should be determined based on
the current value of the instantaneous bias (see [REWEIGHT_BIAS](REWEIGHT_BIAS.md)).

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,3
d3: DISTANCE ATOMS=1,4
d4: DISTANCE ATOMS=2,3
d5: DISTANCE ATOMS=2,4
d6: DISTANCE ATOMS=3,4
rr: RESTRAINT ARG=d1 AT=0.1 KAPPA=10
rbias: REWEIGHT_BIAS TEMP=300

ff: COLLECT_FRAMES ARG=d1,d2,d3,d4,d5,d6 LOGWEIGHTS=rbias STRIDE=5
pca: PCA ARG=ff NLOW_DIM=2 FILE=pca-comp.pdb
DUMPVECTOR ARG=pca,pca_weights FILE=colvar STRIDE=1000
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

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
  keys.add("compulsory","STRIDE","0","the frequency with which to perform this analysis");
  keys.add("optional","FILE","the file on which to output the low dimensional coordinates");
  keys.add("compulsory","FMT","%f","the format to use when outputting the low dimensional coordinates");
  keys.setValueDescription("matrix","the projections of the input coordinates on the PCA components that were found from the covariance matrix");
  keys.needsAction("LOGSUMEXP");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("CONSTANT");
  keys.needsAction("COLLECT");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("CUSTOM");
  keys.needsAction("MATRIX_PRODUCT");
  keys.needsAction("DIAGONALIZE");
  keys.needsAction("VSTACK");
  keys.needsAction("DUMPPDB");
}


PCA::PCA(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Find the argument name
  std::string argn;
  parse("ARG",argn);
  ActionShortcut* as = plumed.getActionSet().getShortcutActionWithLabel( argn );
  if( !as || as->getName()!="COLLECT_FRAMES" ) {
    error("found no COLLECT_FRAMES action with label " + argn );
  }
  // Get the final weights using the logsumexp trick
  readInputLine( getShortcutLabel() + "_weights: LOGSUMEXP ARG=" + argn + "_logweights");
  // Now transpose the collected frames
  readInputLine( getShortcutLabel() + "_dataT: TRANSPOSE ARG=" + argn + "_data");
  // And multiply the transpose by the weights to get the averages
  readInputLine( getShortcutLabel() + "_mean: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_dataT," + getShortcutLabel() + "_weights");
  // Make a matrix of averages
  readInputLine( getShortcutLabel() + "_averages: OUTER_PRODUCT ARG=" + argn + "_ones," + getShortcutLabel() + "_mean");
  // Make a matrix of weights
  ActionWithValue* av2 = plumed.getActionSet().selectWithLabel<ActionWithValue*>( argn + "_data" );
  if( !av2 ) {
    error("count not find data");
  }
  unsigned nones = (av2->copyOutput(0))->getShape()[1];
  std::string ones="1";
  for(unsigned i=1; i<nones; ++i) {
    ones += ",1";
  }
  readInputLine( getShortcutLabel() + "_wones: CONSTANT VALUES=" + ones );
  readInputLine( getShortcutLabel() + "_wmat: OUTER_PRODUCT ARG=" + getShortcutLabel() + "_weights," + getShortcutLabel() + "_wones");
  // And compute the data substract the mean
  readInputLine( getShortcutLabel() + "_diff: CUSTOM ARG=" + argn + "_data," + getShortcutLabel() + "_averages FUNC=(x-y) PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_wdiff: CUSTOM ARG=" + getShortcutLabel() + "_wmat," + getShortcutLabel() + "_diff FUNC=sqrt(x)*y PERIODIC=NO");
  // And the covariance
  readInputLine( getShortcutLabel() + "_wdiffT: TRANSPOSE ARG=" + getShortcutLabel() + "_wdiff");
  readInputLine( getShortcutLabel() + "_covar: MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_wdiffT," + getShortcutLabel() + "_wdiff");
  // Read the dimensionality of the low dimensional space
  unsigned ndim;
  parse("NLOW_DIM",ndim);
  std::string vecstr="1";
  if( ndim<=0 || ndim>nones ) {
    error("cannot generate projection in space of dimension higher than input coordinates");
  }
  for(unsigned i=1; i<ndim; ++i) {
    std::string num;
    Tools::convert( i+1, num );
    vecstr += "," + num;
  }
  readInputLine( getShortcutLabel() + "_eig: DIAGONALIZE ARG=" + getShortcutLabel() + "_covar VECTORS=" + vecstr );
  // Now create a matrix to hold the output data
  std::string outd = "ARG=" + getShortcutLabel() + "_mean";
  for(unsigned i=0; i<ndim; ++i) {
    std::string num;
    Tools::convert( i+1, num );
    outd += "," + getShortcutLabel() + "_eig.vecs-" + num;
  }
  readInputLine( getShortcutLabel() + "_pcaT: VSTACK " + outd );
  readInputLine( getShortcutLabel() + "_pca: TRANSPOSE ARG=" + getShortcutLabel() + "_pcaT");
  // And output it all
  std::string filename, pstride;
  parse("STRIDE",pstride);
  parse("FILE",filename);
  if( filename.length()>0 && av2->getName()=="VSTACK" ) {
    std::vector<std::string> argnames;
    av2->getMatrixColumnTitles( argnames );
    std::string argname_str=argnames[0];
    for(unsigned i=1; i<argnames.size(); ++i) {
      argname_str += "," + argnames[i];
    }
    std::string fmt;
    parse("FMT",fmt);
    if( fmt.length()>0 ) {
      fmt=" FMT=" + fmt;
    }
    readInputLine("DUMPPDB DESCRIPTION=PCA ARG_NAMES=" + argname_str + " ARG=" + getShortcutLabel() + "_pca FILE=" + filename + " STRIDE=" + pstride + fmt );
  } else {
    if( av2->getName()!="COLLECT" ) {
      error("input data should be VSTACK if list of arguments of COLLECT if atom positions");
    }
    ActionAtomistic* rmsdact = plumed.getActionSet().selectWithLabel<ActionAtomistic*>( argn + "_getposx" );
    if( !rmsdact ) {
      error("could not find action that gets positions from trajectory for RMSD");
    }
    std::vector<AtomNumber> atoms( rmsdact->getAbsoluteIndexes() );
    std::string indices;
    Tools::convert( atoms[0].serial(), indices );
    for(unsigned i=1; i<atoms.size(); ++i) {
      std::string jnum;
      Tools::convert( atoms[i].serial(), jnum );
      indices += "," + jnum;
    }
    readInputLine("DUMPPDB DESCRIPTION=PCA ATOM_INDICES=" + indices + " ATOMS=" + getShortcutLabel() + "_pca FILE=" + filename + " STRIDE=" + pstride );
  }
  outd = "ARG=" + getShortcutLabel() + "_eig.vecs-1";
  for(unsigned i=1; i<ndim; ++i) {
    std::string num;
    Tools::convert( i+1, num );
    outd += "," + getShortcutLabel() + "_eig.vecs-" + num;
  }
  readInputLine( getShortcutLabel() + "_eigv: VSTACK " + outd );
  readInputLine( getShortcutLabel() + ": MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_diff," + getShortcutLabel() + "_eigv");
}

}
}
