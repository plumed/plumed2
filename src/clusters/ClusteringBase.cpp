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

namespace PLMD {
namespace clusters {

void ClusteringBase::registerKeywords( Keywords& keys ) {
  matrixtools::MatrixOperationBase::registerKeywords( keys );
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs ");
  keys.setValueDescription("vector","vector with length that is equal to the number of rows in the input matrix.  Elements of this vector are equal to the cluster that each node is a part of");
  keys.addDOI("10.1021/acs.jctc.6b01073");
}

ClusteringBase::ClusteringBase(const ActionOptions&ao):
  Action(ao),
  matrixtools::MatrixOperationBase(ao),
  number_of_cluster(-1) {
  // Do some checks on the input
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) {
    error("input matrix should be square");
  }

  // Now create a value - this holds the data on which cluster each guy is in
  std::vector<std::size_t> shape(1);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  // Build the store here to make sure that next action has all data
  addValue( shape );
  setNotPeriodic();
  // Resize local variables
  which_cluster.resize( getPntrToArgument(0)->getShape()[0] );
  cluster_sizes.resize( getPntrToArgument(0)->getShape()[0] );
  log<<"  Bibliography "<<plumed.cite("10.1021/acs.jctc.6b01073")<<"\n";
}

void ClusteringBase::retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list ) {
  // Make sure we have the edges stored
  std::vector<std::pair<unsigned,unsigned> > pairs;
  std::vector<double> vals;
  unsigned nedge;
  getPntrToArgument(0)->retrieveEdgeList( nedge, pairs, vals );
  // Currently everything has zero neighbors
  for(unsigned i=0; i<nneigh.size(); ++i) {
    nneigh[i]=0;
  }
  // Resize the adjacency list if it is needed
  if( adj_list.ncols()!=getPntrToArgument(0)->getNumberOfColumns() ) {
    unsigned nrows = getPntrToArgument(0)->getShape()[0];
    adj_list.resize( nrows, getPntrToArgument(0)->getNumberOfColumns() );
  }

  // And set up the adjacency list
  for(unsigned i=0; i<nedge; ++i) {
    // Store if atoms are connected
    unsigned j=pairs[i].first, k=pairs[i].second;
    if( j==k ) {
      continue;
    }
    adj_list(j,nneigh[j])=k;
    nneigh[j]++;
    adj_list(k,nneigh[k])=j;
    nneigh[k]++;
  }
}

void ClusteringBase::calculate() {
  // All the clusters have zero size initially
  for(unsigned i=0; i<cluster_sizes.size(); ++i) {
    cluster_sizes[i].first=0;
    cluster_sizes[i].second=i;
  }
  // Do the clustering bit
  performClustering();
  // Order the clusters in the system by size (this returns ascending order )
  std::sort( cluster_sizes.begin(), cluster_sizes.end() );
  // Set the elements of the value to the cluster identies
  for(unsigned i=0; i<cluster_sizes.size(); ++i) {
    double this_size = static_cast<double>(cluster_sizes.size()-i);
    for(unsigned j=0; j<cluster_sizes.size(); ++j) {
      if( which_cluster[j]==cluster_sizes[i].second ) {
        getPntrToValue()->set( j, this_size );
      }
    }
  }
}

void ClusteringBase::apply() {
  if( getPntrToComponent(0)->forcesWereAdded() ) {
    error("forces on clustering actions cannot work as clustering is not differentiable");
  }
}

}
}
