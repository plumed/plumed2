/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
namespace adjmat {

void ClusteringBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.add("compulsory","MATRIX","the label of the adjacency matrix we would like to analyse");
}

ClusteringBase::ClusteringBase(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao),
  number_of_cluster(-1)
{
  std::vector<Value*> mat; parseArgumentList("MATRIX",mat);
  if( mat.size()!=1 ) error("should pass only one matrix to clustering base");
  // Check if we have a matrix
  if( mat[0]->getRank()!=2 ) error("input should be a matrix");
  // Check if matrix has same number of rows and columns
  if( mat[0]->getShape()[0]!=mat[0]->getShape()[1] ) error("matrix should be square");
  log.printf("  finding clusters in matrix labelled %s \n", mat[0]->getName().c_str() );
  // Request the argument
  requestArguments(mat, false); checkRead();
  // Now create a value - this holds the data on which cluster each guy is in
  std::vector<unsigned> shape(1); shape[0]=mat[0]->getShape()[0];
  // Build the store here to make sure that next action has all data
  addValue( shape ); setNotPeriodic(); getPntrToOutput(0)->alwaysStoreValues();
  // Resize local variables
  which_cluster.resize( mat[0]->getShape()[0] ); cluster_sizes.resize( mat[0]->getShape()[0] );
  // Create a group for this action
  const auto m=plumed.getAtoms().getAllGroups().find(mat[0]->getPntrToAction()->getLabel());
  plumed.getAtoms().insertGroup( getLabel(), m->second );
}

void ClusteringBase::retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list ) {
  // Currently everything has zero neighbors
  for(unsigned i=0; i<nneigh.size(); ++i) nneigh[i]=0;

  // And set up the adjacency list
  Value* mat = getPntrToArgument(0); unsigned nrows = mat->getShape()[0];
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
  nedge=0; std::vector<double> myvals( getNumberOfComponents() );
  if( getPntrToArgument(0)->getNumberOfValues(getLabel())>edge_list.size() ) error("adjacency lists are not large enough, increase maxconnections");

  Value* mat = getPntrToArgument(0); unsigned nrows = mat->getShape()[0];
  for(unsigned i=0; i<mat->getNumberOfValues(getLabel()); ++i) {
    // Check if atoms are connected
    if( mat->get(i)<epsilon ) continue ;

    unsigned j = std::floor( i / nrows ); unsigned k = i%nrows;
    if( k>j ) continue ;  // Ensures each connection is only stored once

    edge_list[nedge].first = j; edge_list[nedge].second = k; nedge++;
  }
}

void ClusteringBase::calculate() {
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

}
}
