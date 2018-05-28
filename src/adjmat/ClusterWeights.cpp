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
#include "core/Atoms.h"
#include "multicolvar/MultiColvarBase.h"
#include "ClusteringBase.h"

//+PLUMEDOC CONCOMP CLUSTER_PROPERTIES
/*
Calculate properties of the distribution of some quantities that are part of a connected component

This collective variable was developed for looking at nucleation phenomena, where you are
interested in using studying the behavior of atoms in small aggregates or nuclei.  In these sorts of
problems you might be interested in the degree the atoms in a nucleus have adopted their crystalline
structure or (in the case of heterogenous nucleation of a solute from a solvent) you might be
interested in how many atoms are present in the largest cluster \cite tribello-clustering.

\par Examples

The input below calculates the coordination numbers of atoms 1-100 and then computes the an adjacency
matrix whose elements measures whether atoms \f$i\f$ and \f$j\f$ are within 0.55 nm of each other.  The action
labelled dfs then treats the elements of this matrix as zero or ones and thus thinks of the matrix as defining
a graph.  This dfs action then finds the largest connected component in this graph.  The sum of the coordination
numbers for the atoms in this largest connected component are then computed and this quantity is output to a colvar
file.  The way this input can be used is described in detail in \cite tribello-clustering.

\plumedfile
lq: COORDINATIONNUMBER SPECIES=1-100 SWITCH={CUBIC D_0=0.45  D_MAX=0.55} LOWMEM
cm: CONTACT_MATRIX ATOMS=lq  SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
dfs: DFSCLUSTERING MATRIX=cm
clust1: CLUSTER_PROPERTIES CLUSTERS=dfs CLUSTER=1 SUM
PRINT ARG=clust1.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class ClusterWeights :
  public ActionWithArguments,
  public ActionWithValue {
private:
/// The cluster we are looking for
  unsigned clustr;
/// The forces
  std::vector<double> forcesToApply;
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ClusterWeights(const ActionOptions&);
/// The number of derivatives
  unsigned getNumberOfDerivatives() const ;
/// Work out what needs to be done in this action
  void buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags );
/// Do the calculation
  void calculate();
/// We can use ActionWithVessel to run all the calculation
  void performTask( const unsigned&, MultiValue& ) const ;
  void apply() {}
};

PLUMED_REGISTER_ACTION(ClusterWeights,"CLUSTER_WEIGHTS")
PLUMED_REGISTER_SHORTCUT(ClusterWeights,"CLUSTER_WEIGHTS")
PLUMED_REGISTER_SHORTCUT(ClusterWeights,"CLUSTER_PROPERTIES")
PLUMED_REGISTER_SHORTCUT(ClusterWeights,"CLUSTER_NATOMS")
PLUMED_REGISTER_SHORTCUT(ClusterWeights,"CLUSTER_WITHSURFACE")
PLUMED_REGISTER_SHORTCUT(ClusterWeights,"CLUSTER_DIAMETER")

void ClusterWeights::shortcutKeywords( Keywords& keys ) {
  keys.add("optional","ARG","calculate the sum of the arguments calculated by this action for the cluster");
  keys.add("optional","RCUT_SURF","");
  multicolvar::MultiColvarBase::shortcutKeywords( keys );
}

void ClusterWeights::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                     const std::map<std::string,std::string>& keys,
                                     std::vector<std::vector<std::string> >& actions ) {
  if( words[0]!="CLUSTER_DIAMETER" ) {
    std::vector<std::string> weights;
    if( words[0]=="CLUSTER_NATOMS" ) weights.push_back( lab + "_weights:" );
    else if( words[0]=="CLUSTER_WITHSURFACE" ) weights.push_back( lab + "_wnosurf:");
    else weights.push_back( lab + ":" );
    weights.push_back("CLUSTER_WEIGHTS");
    for(unsigned i=1; i<words.size(); ++i) weights.push_back( words[i] );
    actions.push_back( weights );
  }
  if( words[0]=="CLUSTER_PROPERTIES" ) {
    multicolvar::MultiColvarBase::expandFunctions( lab, keys.find("ARG")->second, lab, words, keys, actions );
  } else if( words[0]=="CLUSTER_NATOMS" ) {
    std::vector<std::string> nat_str; nat_str.push_back( lab + ":" ); nat_str.push_back("COMBINE");
    nat_str.push_back("ARG=" + lab + "_weights"); nat_str.push_back("PERIODIC=NO"); actions.push_back( nat_str );
  } else if( words[0]=="CLUSTER_WITHSURFACE" ) {
    // Contact matrix
    std::vector<std::string> contact_mat; contact_mat.push_back( lab + "_cmat:"); contact_mat.push_back("CONTACT_MATRIX");
    contact_mat.push_back("GROUP=" + lab + "_wnosurf"); contact_mat.push_back("SWITCH=" + keys.find("RCUT_SURF")->second );
    actions.push_back( contact_mat );
    // Matrix of products of cluster weights
    std::vector<std::string> cweight_mat; cweight_mat.push_back( lab + "_cwmat:"); cweight_mat.push_back("CUSTOM_MATRIX");
    cweight_mat.push_back("GROUP1=" + lab + "_wnosurf"); cweight_mat.push_back("FUNC=max");
    actions.push_back( cweight_mat );
    // Product of matrices
    std::vector<std::string> prod_mat; prod_mat.push_back( lab + "_pmat:"); prod_mat.push_back("MATHEVAL");
    prod_mat.push_back("ARG1=" + lab + "_cmat.w"); prod_mat.push_back("ARG2=" + lab + "_cwmat"); prod_mat.push_back("FUNC=x*y");
    prod_mat.push_back("PERIODIC=NO"); actions.push_back( prod_mat );
    // Clusters
    std::vector<std::string> cluste; cluste.push_back( lab + "_clust:" ); cluste.push_back("DFSCLUSTERING");
    cluste.push_back("MATRIX=" + lab + "_pmat"); actions.push_back( cluste );
    // Final weights
    std::vector<std::string> fweights; fweights.push_back( lab + ":" ); fweights.push_back("CLUSTER_WEIGHTS");
    fweights.push_back("CLUSTERS=" + lab + "_clust"); fweights.push_back("CLUSTER=1");
    actions.push_back( fweights );
  } else if( words[0]=="CLUSTER_DIAMETER" ) {
    // Distance matrix
    std::vector<std::string> distance_mat; distance_mat.push_back( lab + "_dmat:" ); distance_mat.push_back("DISTANCE_MATRIX");
    distance_mat.push_back("GROUP=" + keys.find("ARG")->second ); actions.push_back( distance_mat );
    // Matrix of bonds in cluster
    std::vector<std::string> bonds_mat; bonds_mat.push_back( lab + "_bmat:"); bonds_mat.push_back("DOTPRODUCT_MATRIX");
    bonds_mat.push_back("GROUP1=" + keys.find("ARG")->second ); actions.push_back( bonds_mat );
    // Product of matrices
    std::vector<std::string> dcls_mat; dcls_mat.push_back( lab + "_dcls:"); dcls_mat.push_back("MATHEVAL");
    dcls_mat.push_back("ARG1=" + lab + "_dmat.w" ); dcls_mat.push_back("ARG2=" + lab + "_bmat");
    dcls_mat.push_back("FUNC=x*y"); dcls_mat.push_back("PERIODIC=NO"); actions.push_back( dcls_mat );
    // And take the highest value
    std::vector<std::string> maxrad; maxrad.push_back( lab + ":"); maxrad.push_back("HIGHEST");
    maxrad.push_back("ARG=" + lab + "_dcls"); actions.push_back( maxrad );
  }
}

void ClusterWeights::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.remove("ARG");
  ActionWithValue::registerKeywords( keys ); keys.remove("NUMERICAL_DERIVATIVES");
  keys.add("compulsory","CLUSTERS","the label of the action that does the clustering");
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
}

ClusterWeights::ClusterWeights(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao)
{
  // Read in the clustering object
  std::vector<Value*> clusters; parseArgumentList("CLUSTERS",clusters);
  if( clusters.size()!=1 ) error("should pass only one matrix to clustering base");
  ClusteringBase* cc = dynamic_cast<ClusteringBase*>( clusters[0]->getPntrToAction() );
  if( !cc ) error("input to CLUSTERS keyword should be a clustering action");
  // Request the arguments
  requestArguments( clusters, false );
  // Now create the value
  std::vector<unsigned> shape(1); shape[0]=clusters[0]->getShape()[0];
  addValue( shape ); setNotPeriodic(); getPntrToOutput(0)->alwaysStoreValues();
  // And the tasks
  for(unsigned i=0; i<shape[0]; ++i) addTaskToList(i);
  // Find out which cluster we want
  parse("CLUSTER",clustr);
  if( clustr<1 ) error("cannot look for a cluster larger than the largest cluster");
  if( clustr>clusters[0]->getShape()[0] ) error("cluster selected is invalid - too few atoms in system");
  log.printf("  atoms in %dth largest cluster calculated by %s are equal to one \n",clustr, cc->getLabel().c_str() );
  // Create a group for this action so we can associate atoms to these weights easily
  const auto m=plumed.getAtoms().getAllGroups().find(clusters[0]->getPntrToAction()->getLabel());
  plumed.getAtoms().insertGroup( getLabel(), m->second );
}

void ClusterWeights::buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
  plumed_assert( getPntrToArgument(0)->valueHasBeenSet() ); actionsThatSelectTasks.push_back( getLabel() );
  for(unsigned i=0; i<getPntrToArgument(0)->getShape()[0]; ++i) {
    if( fabs(getPntrToArgument(0)->get(i)-clustr)<epsilon ) tflags[i]=1;
  }
}

unsigned ClusterWeights::getNumberOfDerivatives() const {
  return 0;
}

void ClusterWeights::calculate() {
  runAllTasks();
}

void ClusterWeights::performTask( const unsigned& current, MultiValue& myvals ) const {
  myvals.addValue( getPntrToOutput(0)->getPositionInStream(), 1.0 );
}

}
}
