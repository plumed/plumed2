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
#include "ClusteringBase.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace clusters {

void ClusteringBase::registerKeywords( Keywords& keys ) {
  matrixtools::ActionWithInputMatrices::registerKeywords( keys ); keys.use("ARG");
}

ClusteringBase::ClusteringBase(const ActionOptions&ao):
  Action(ao),
  matrixtools::ActionWithInputMatrices(ao),
  number_of_cluster(-1)
{
  // Do some checks on the input
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) error("input matrix should be square"); 

  // Now create a value - this holds the data on which cluster each guy is in
  std::vector<unsigned> shape(1); shape[0]=getPntrToArgument(0)->getShape()[0];
  // Build the store here to make sure that next action has all data
  addValue( shape ); setNotPeriodic(); getPntrToOutput(0)->alwaysStoreValues();
  // Resize local variables
  which_cluster.resize( getPntrToArgument(0)->getShape()[0] ); cluster_sizes.resize( getPntrToArgument(0)->getShape()[0] );
  // Create a group for this action
  const auto m=plumed.getAtoms().getAllGroups().find(getPntrToArgument(0)->getPntrToAction()->getLabel());
  plumed.getAtoms().insertGroup( getLabel(), m->second );
}

void ClusteringBase::retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list ) {
  Value* mat = getPntrToArgument(0); unsigned nrows = mat->getShape()[0];
  // Currently everything has zero neighbors
  for(unsigned i=0; i<nneigh.size(); ++i) nneigh[i]=0;
  // Resize the adjacency list if it is needed
  if( adj_list.ncols()!=mat->getNumberOfColumns() ) adj_list.resize( nrows, mat->getNumberOfColumns() );

  // And set up the adjacency list
  for(unsigned i=0; i<mat->getNumberOfValues(getLabel()); ++i) {
    // Check if atoms are connected
    if( mat->get(i)<epsilon ) continue ;

    unsigned j = std::floor( i / nrows ); unsigned k = i%nrows;
    if( k>j ) continue ;  // Ensures each connection is only stored once

    if( nneigh[j]>=adj_list.ncols() || nneigh[k]>=adj_list.ncols() ) error("adjacency lists are not large enough, increase maxconnections");
    // Store if atoms are connected
    adj_list(k,nneigh[k])=j; nneigh[k]++;
    adj_list(j,nneigh[j])=k; nneigh[j]++;
  }
}

void ClusteringBase::retrieveEdgeList( unsigned& nedge, std::vector<std::pair<unsigned,unsigned> >& edge_list ) {
  nedge=0; Value* mat = getPntrToArgument(0); unsigned nrows = mat->getShape()[0];
  if( edge_list.size()!=nrows*mat->getNumberOfColumns() ) edge_list.resize( nrows*mat->getNumberOfColumns() );

  for(unsigned i=0; i<mat->getNumberOfValues(getLabel()); ++i) {
    // Check if atoms are connected
    if( mat->get(i)<epsilon ) continue ;

    unsigned j = std::floor( i / nrows ); unsigned k = i%nrows;
    if( k>j ) continue ;  // Ensures each connection is only stored once

    edge_list[nedge].first = j; edge_list[nedge].second = k; nedge++;
  }
}

void ClusteringBase::completeMatrixOperations() {
  // All the clusters have zero size initially
  for(unsigned i=0; i<cluster_sizes.size(); ++i) { cluster_sizes[i].first=0; cluster_sizes[i].second=i; }
  // Do the clustering bit
  performClustering();
  // Order the clusters in the system by size (this returns ascending order )
  std::sort( cluster_sizes.begin(), cluster_sizes.end() );
  // Set the elements of the value to the cluster identies
  for(unsigned i=0; i<cluster_sizes.size(); ++i) {
    double this_size = static_cast<double>(cluster_sizes.size()-i);
    for(unsigned j=0; j<cluster_sizes.size(); ++j) {
      if( which_cluster[j]==cluster_sizes[i].second ) getPntrToValue()->set( j, this_size );
    }
  }
}

void ClusteringBase::apply() {
  if( getPntrToOutput(0)->forcesWereAdded() ) error("forces on clustering actions cannot work as clustering is not differentiable"); 
}

}
}
