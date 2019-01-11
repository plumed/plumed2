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
#include "AdjacencyMatrixBase.h"
#include "AdjacencyMatrixVessel.h"

namespace PLMD {
namespace adjmat {

void ClusteringBase::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrix::registerKeywords( keys );
}

ClusteringBase::ClusteringBase(const ActionOptions&ao):
  Action(ao),
  ActionWithInputMatrix(ao),
  number_of_cluster(-1)
{
  if( getAdjacencyVessel() ) {
    cluster_sizes.resize(getNumberOfNodes()); which_cluster.resize(getNumberOfNodes());
    if( getNumberOfNodeTypes()!=1 ) error("should only be running clustering with one base multicolvar in function");
    if( !getAdjacencyVessel()->undirectedGraph() ) error("input contact matrix is incompatible with clustering");
  }
  if( keywords.exists("MATRIX") ) {
    std::vector<AtomNumber> fake_atoms; setupMultiColvarBase( fake_atoms );
  }
}

void ClusteringBase::turnOnDerivatives() {
  // Check base multicolvar isn't density probably other things shouldn't be allowed here as well
  if( (getAdjacencyVessel()->getMatrixAction())->getNumberOfBaseMultiColvars()>0 ) {
    if( getBaseMultiColvar(0)->isDensity() ) error("DFS clustering cannot be differentiated if base multicolvar is DENSITY");
  }

  // Ensure that derivatives are turned on in base classes
  ActionWithInputMatrix::turnOnDerivatives();
}

void ClusteringBase::calculate() {
  // All the clusters have zero size initially
  for(unsigned i=0; i<cluster_sizes.size(); ++i) { cluster_sizes[i].first=0; cluster_sizes[i].second=i; }
  // Do the clustering bit
  performClustering();
  // Order the clusters in the system by size (this returns ascending order )
  std::sort( cluster_sizes.begin(), cluster_sizes.end() );
}

void ClusteringBase::retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const {
  unsigned n=0; myatoms.resize( cluster_sizes[cluster_sizes.size() - clust].first );
  for(unsigned i=0; i<getNumberOfNodes(); ++i) {
    if( which_cluster[i]==cluster_sizes[cluster_sizes.size() - clust].second ) { myatoms[n]=i; n++; }
  }
}

bool ClusteringBase::areConnected( const unsigned& iatom, const unsigned& jatom ) const {
  return getAdjacencyVessel()->nodesAreConnected( iatom, jatom );
}

double ClusteringBase::getCutoffForConnection() const {
  return getAdjacencyVessel()->getCutoffForConnection();
}

}
}
