/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014,2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#ifndef __PLUMED_adjmat_DFSBase_h
#define __PLUMED_adjmat_DFSBase_h

#include "ActionWithInputMatrix.h"
#include "multicolvar/MultiColvarFunction.h"
#include "multicolvar/AtomValuePack.h"

namespace PLMD {
namespace adjmat {

class DFSBase : public ActionWithInputMatrix {
private:
/// Vector that stores the sizes of the current set of clusters
  std::vector< std::pair<unsigned,unsigned> > cluster_sizes;
/// Vector that identifies the cluster each atom belongs to
  std::vector<unsigned> which_cluster;
/// Used to identify the cluster we are working on 
  int number_of_cluster;
#ifdef __PLUMED_HAS_BOOST
/// The list of edges in the graph
  std::vector<std::pair<unsigned,unsigned> > edge_list;
#else 
/// The number of neighbors each atom has
  std::vector<unsigned> nneigh;
/// The adjacency list
  Matrix<unsigned> adj_list;
/// The color that tells us whether a node has been visited
  std::vector<unsigned> color;
/// The recursive function at the heart of this method
  int explore( const unsigned& index );
#endif
protected:
/// Get the number of clusters that have been found
  unsigned getNumberOfClusters() const ;
/// Get the atoms in one of the clusters
  void retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const ;
/// Do the clustering of the dat
  void performClustering();
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DFSBase(const ActionOptions&);
/// Required as we have to be able to deal with vectors
  unsigned getNumberOfQuantities();
/// This checks whether derivatives can be computed given the base multicolvar
  void turnOnDerivatives();
};

inline
unsigned DFSBase::getNumberOfClusters() const {
  return number_of_cluster + 1;
}

}
}

#endif
