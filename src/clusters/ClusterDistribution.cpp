/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "ClusteringBase.h"

//+PLUMEDOC CONCOMP CLUSTER_DISTRIBUTION
/*
Calculate functions of the distribution of properties in your connected components.

This action allows you to calculate the number of atoms in each of the connected components that
were detected when you performed [DFSCLUSTERING](DFSCLUSTERING.md) on one of the adjacency matrices
computed using the [admat](module_adjmat.md).  The following example illustrates how you can use
this to compute the number of atoms in each of the identified clusters:

```plumed
cm: CONTACT_MATRIX GROUP=1-100 SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
dfs: DFSCLUSTERING ARG=cm
clust: CLUSTER_DISTRIBUTION CLUSTERS=dfs
PRINT ARG=clust FILE=colvar STRIDE=1
```

The output from the CLUSTER_DISTRIBUTION action here is a vector with 100 elements.  The first element of this vector
is the number of atoms in the largest connected component, the second element of the vector is the number of atoms in the
second largest connected component and so on (many elements of the vector will be zero).

As illustrated in the inputs below there are multiple shortcuts that allow you to probe the distribution of cluster sizes.
For example, the following input calculates how many clusters with more than 10 atoms are present.

```plumed
cm: CONTACT_MATRIX GROUP=1-100 SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
dfs: DFSCLUSTERING ARG=cm
clust: CLUSTER_DISTRIBUTION CLUSTERS=dfs MORE_THAN={RATIONAL D_0=10 R_0=0.0001}
PRINT ARG=clust.morethan FILE=colvar STRIDE=1
```

By a similar logic you can compute the number of clusters that were identified as follows:

```plumed
cm: CONTACT_MATRIX GROUP=1-100 SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
dfs: DFSCLUSTERING ARG=cm
clust: CLUSTER_DISTRIBUTION CLUSTERS=dfs LESS_THAN={RATIONAL R_0=0.0001}
nc: CUSTOM ARG=clust.lessthan FUNC=100-x PERIODIC=NO
PRINT ARG=nc FILE=colvar STRIDE=1
```

This input calculates the number of zeros in the vector output by the CLUSTER_DISTRIBUTION action.  If we subtract the number of
non-zero elements from the size of the vector we then get the number of clusters.

Lastly, you can calculate the number of clusters that are within a certain range or a set of ranges using the BETWEEN and HISGTOGRAM keywords
as indicated below:

```plumed
cm: CONTACT_MATRIX GROUP=1-100 SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
dfs: DFSCLUSTERING ARG=cm
clust: CLUSTER_DISTRIBUTION ...
   CLUSTERS=dfs
   BETWEEN={GAUSSIAN LOWER=5 UPPER=6 SMEAR=0.5}
   HISTOGRAM={GAUSSIAN LOWER=6 UPPER=10 NBINS=4 SMEAR=0.5}
...
PRINT ARG=clust.* FILE=colvar STRIDE=1
```

This input will output 5 quantities:

- `clust.between` tells you the number of connected components that contain between 5 and 6 atoms.
- `clust.between-1` tells you the number of connected components that contain between 6 and 7 atoms.
- `clust.between-2` tells you the number of connected components that contain between 7 and 8 atoms.
- `clust.between-3` tells you the number of connected components that contain between 8 and 9 atoms.
- `clust.between-4` tells you the number of connected components that contain between 9 and 10 atoms.

If you expand the inputs above you can find more details on how these quantities are calculated.

The input provided below shows just how this can be used to perform quite complicated calculations.  The input
calculates the local q6 Steinhardt parameter on each atom.  The coordination number
that atoms with a high value for the local q6 Steinhardt parameter have with other atoms that have a high
value for the local q6 Steinhardt parameter is then computed.  A contact matrix is then computed that measures
whether atoms atoms $i$ and $j$ have a high value for this coordination number and if they are within
3.6 nm of each other.  The connected components of this matrix are then found using a depth first clustering
algorithm on the corresponding graph. The number of components in this graph that contain more than 27 atoms is then computed.
An input similar to this one was used to analyze the formation of a polycrystal of GeTe from amorphous GeTe in the paper cited below

```plumed
q6: Q6 SPECIES=1-300 SWITCH={GAUSSIAN D_0=5.29 R_0=0.01 D_MAX=5.3}
lq6: LOCAL_Q6 SPECIES=q6 SWITCH={GAUSSIAN D_0=5.29 R_0=0.01 D_MAX=5.3}
flq6: MORE_THAN ARG=lq6 SWITCH={GAUSSIAN D_0=0.19 R_0=0.01 D_MAX=0.2}
cc: COORDINATIONNUMBER SPECIES=1-300 SWITCH={GAUSSIAN D_0=3.59 R_0=0.01 D_MAX=3.6}
fcc: MORE_THAN ARG=cc SWITCH={GAUSSIAN D_0=5.99 R_0=0.01 D_MAX=6.0}
mat: CONTACT_MATRIX GROUP=1-300 SWITCH={GAUSSIAN D_0=3.59 R_0=0.01 D_MAX=3.6}
dfs: DFSCLUSTERING ARG=mat
nclust: CLUSTER_DISTRIBUTION ...
   CLUSTERS=dfs WEIGHTS=fcc
   MORE_THAN={GAUSSIAN D_0=26.99 R_0=0.01 D_MAX=27}
...
PRINT ARG=nclust.* FILE=colvar
```

Notice how the WEIGHTS keyword is used here in the input to CLUSTER_DISTRIBUTION. By using this keyword here we ensure that the size of each
cluster is the sum of the components of the vector `fcc` that are part of the connected component.  The `size` of each cluster that is output
is thus the number of atoms in that cluster that have a COORDINATIONNUMBER that is greater than 4.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace clusters {

class ClusterDistribution :
  public ActionWithArguments,
  public ActionWithValue {
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ClusterDistribution(const ActionOptions&);
/// The number of derivatives
  unsigned getNumberOfDerivatives() override ;
/// Do the calculation
  void calculate() override;
/// We can use ActionWithVessel to run all the calculation
  void apply() override {}
};

PLUMED_REGISTER_ACTION(ClusterDistribution,"CLUSTER_DISTRIBUTION_CALC")

void ClusterDistribution::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.setDisplayName("CLUSTER_DISTRIBUTION");
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addInputKeyword("compulsory","CLUSTERS","vector","the label of the action that does the clustering");
  keys.addInputKeyword("optional","WEIGHTS","vector","use the vector of values calculated by this action as weights rather than giving each atom a unit weight");
  keys.setValueDescription("vector","a vector containing the sum of a atomic-cv that is calculated for each of the identified clusters");
  keys.addDOI("https://doi.org/10.1021/acs.jctc.6b01073");
}

ClusterDistribution::ClusterDistribution(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao) {
  // Read in the clustering object
  std::vector<Value*> clusters;
  parseArgumentList("CLUSTERS",clusters);
  if( clusters.size()!=1 ) {
    error("should pass only one matrix to clustering base");
  }
  ClusteringBase* cc = dynamic_cast<ClusteringBase*>( clusters[0]->getPntrToAction() );
  if( !cc ) {
    error("input to CLUSTERS keyword should be a clustering action");
  }
  std::vector<Value*> weights;
  parseArgumentList("WEIGHTS",weights);
  if( weights.size()==0 ) {
    log.printf("  using unit weights in cluster distribution \n");
  } else if( weights.size()==1 ) {
    if( weights[0]->getRank()!=1 ) {
      error("input weights has wrong shape");
    }
    if( weights[0]->getShape()[0]!=clusters[0]->getShape()[0] ) {
      error("mismatch between number of weights and number of atoms");
    }
    log.printf("  using weights from action with label %s in cluster distribution \n", weights[0]->getName().c_str() );
    clusters.push_back( weights[0] );
  } else {
    error("should have only one argument for weights \n");
  }
  // Request the arguments
  requestArguments( clusters );
  // Now create the value
  std::vector<std::size_t> shape(1);
  shape[0]=clusters[0]->getShape()[0];
  addValue( shape );
  setNotPeriodic();
}

unsigned ClusterDistribution::getNumberOfDerivatives() {
  return 0;
}

void ClusterDistribution::calculate() {
  plumed_assert( getPntrToArgument(0)->valueHasBeenSet() );
  if( getNumberOfArguments()>1 ) {
    plumed_assert( getPntrToArgument(1)->valueHasBeenSet() );
  }
  double csize = getPntrToArgument(0)->get(0);
  for(unsigned i=1; i<getPntrToArgument(0)->getShape()[0]; ++i) {
    if( getPntrToArgument(0)->get(i)>csize ) {
      csize = getPntrToArgument(0)->get(i);
    }
  }
  unsigned ntasks = static_cast<unsigned>( csize );
  for(unsigned i=0; i<ntasks; ++i) {
    for(unsigned j=0; j<getPntrToArgument(0)->getShape()[0]; ++j) {
      if( fabs(getPntrToArgument(0)->get(j)-i)<epsilon ) {
        if( getNumberOfArguments()==2 ) {
          getPntrToValue()->add( i, getPntrToArgument(1)->get(j) );
        } else {
          getPntrToValue()->add( i, 1.0 );
        }
      }
    }
  }
}

class ClusterDistributionShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  ClusterDistributionShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ClusterDistributionShortcut,"CLUSTER_DISTRIBUTION")

void ClusterDistributionShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","CLUSTERS","the label of the action that does the clustering");
  keys.add("optional","WEIGHTS","use the vector of values calculated by this action as weights rather than giving each atom a unit weight");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.reset_style("MIN","hidden");
  keys.reset_style("MAX","hidden");
  keys.reset_style("ALT_MIN","hidden");
  keys.reset_style("HIGHEST","hidden");
  keys.reset_style("LOWEST","hidden");
  keys.reset_style("MEAN","hidden");
  keys.reset_style("SUM","hidden");
  keys.addActionNameSuffix("_CALC");
}

ClusterDistributionShortcut::ClusterDistributionShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  readInputLine( getShortcutLabel() + ": CLUSTER_DISTRIBUTION_CALC " + convertInputLineToString() );
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(),  getShortcutLabel(),  "", keymap, this );
}

}
}
