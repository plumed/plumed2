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
#ifndef __PLUMED_adjmat_ClusteringBase_h
#define __PLUMED_adjmat_ClusteringBase_h

#include "ActionWithInputMatrix.h"
#include "multicolvar/AtomValuePack.h"

namespace PLMD {
namespace adjmat {

class ClusteringBase : public ActionWithInputMatrix {
protected:
/// Vector that stores the sizes of the current set of clusters
  std::vector< std::pair<unsigned,unsigned> > cluster_sizes;
/// Used to identify the cluster we are working on
  int number_of_cluster;
/// Vector that identifies the cluster each atom belongs to
  std::vector<unsigned> which_cluster;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ClusteringBase(const ActionOptions&);
/// This checks whether derivatives can be computed given the base multicolvar
  void turnOnDerivatives();
/// Are these two atoms connected
  bool areConnected( const unsigned& iatom, const unsigned& jatom ) const ;
/// Do the calculation
  void calculate();
/// Do the clustering
  virtual void performClustering()=0;
/// Get the number of clusters that have been found
  unsigned getNumberOfClusters() const ;
/// Get the atoms in one of the clusters
  virtual void retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const ;
/// Do nothing for apply here
  void apply() {}
/// Get the cutoff
  virtual double getCutoffForConnection() const ;
};

inline
unsigned ClusteringBase::getNumberOfClusters() const {
  return number_of_cluster + 1;
}

}
}

#endif
