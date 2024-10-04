/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "ClusteringBase.h"

//+PLUMEDOC CONCOMP CLUSTER_WEIGHTS
/*
Setup a vector that has one for all the atoms that form part of the cluster of interest and that has zero for all other atoms.

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace clusters {

class ClusterWeights :
  public ActionWithArguments,
  public ActionWithValue {
private:
/// The cluster we are looking for
  unsigned clustr;
/// The forces
  std::vector<double> forcesToApply;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ClusterWeights(const ActionOptions&);
/// The number of derivatives
  unsigned getNumberOfDerivatives() override ;
/// Do the calculation
  void calculate() override ;
///
  void apply() override {}
};

PLUMED_REGISTER_ACTION(ClusterWeights,"CLUSTER_WEIGHTS")

void ClusterWeights::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.remove("ARG");
  ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.add("compulsory","CLUSTERS","the label of the action that does the clustering");
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
  keys.addFlag("LOWMEM",false,"this flag does nothing and is present only to ensure back-compatibility");
  // keys.add("hidden","FROM_PROPERTIES","indicates that this is created from CLUSTER_PROPERTIES shortcut");
  keys.setValueDescription("vector with elements that are one if the atom of interest is part of the required cluster and zero otherwise");
}

ClusterWeights::ClusterWeights(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao) {
  bool lowmem;
  parseFlag("LOWMEM",lowmem);
  if( lowmem ) {
    warning("LOWMEM flag is deprecated and is no longer required for this action");
  }
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
  // Request the arguments
  requestArguments( clusters );
  // Now create the value
  std::vector<unsigned> shape(1);
  shape[0]=clusters[0]->getShape()[0];
  addValue( shape );
  setNotPeriodic();
  getPntrToComponent(0)->buildDataStore();
  // Find out which cluster we want
  parse("CLUSTER",clustr);
  if( clustr<1 ) {
    error("cannot look for a cluster larger than the largest cluster");
  }
  if( clustr>clusters[0]->getShape()[0] ) {
    error("cluster selected is invalid - too few atoms in system");
  }
  log.printf("  atoms in %dth largest cluster calculated by %s are equal to one \n",clustr, cc->getLabel().c_str() );
}

unsigned ClusterWeights::getNumberOfDerivatives() {
  return 0;
}

void ClusterWeights::calculate() {
  plumed_assert( getPntrToArgument(0)->valueHasBeenSet() );
  for(unsigned i=0; i<getPntrToArgument(0)->getShape()[0]; ++i) {
    if( fabs(getPntrToArgument(0)->get(i)-clustr)<epsilon ) {
      getPntrToComponent(0)->set( i, 1.0 );
    }
  }
}

}
}
