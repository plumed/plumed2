/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#ifndef __PLUMED_adjmat_AdjacencyMatrixVessel_h
#define __PLUMED_adjmat_AdjacencyMatrixVessel_h

#include "vesselbase/StoreDataVessel.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace adjmat {

class AdjacencyMatrixBase;

// One school of thought would have it that it makes more sense to
// have the functionality contained within this class in AdjacencyMatrixBase
// I have not done this as I can inherit many useful things from StoreDataVessel
// If I put this functionality within AdjacencyMatrixBase I would have to reimplement
// these features.

class AdjacencyMatrixVessel : public vesselbase::StoreDataVessel {
  friend class AdjacencyMatrixBase;
  friend class ActionWithInputMatrix;
private:
/// Pointer to underlying action
  AdjacencyMatrixBase* function;
/// Is the matrix symmetric and are we calculating hbonds
  bool symmetric, hbonds;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit AdjacencyMatrixVessel( const vesselbase::VesselOptions& );
/// Get the underlying adjacency matrix action object
  AdjacencyMatrixBase* getMatrixAction();
/// Is an element of the matrix currently active
  bool matrixElementIsActive( const unsigned& ielem, const unsigned& jelem ) const ;
/// Get the index that a particular element is stored in from the matrix indices
  unsigned getStoreIndexFromMatrixIndices( const unsigned& ielem, const unsigned& jelem ) const ;
/// Get the adjacency matrix
  void retrieveMatrix( DynamicList<unsigned>& myactive_elements, Matrix<double>& mymatrix );
/// Get the neighbour list based on the adjacency matrix
  void retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list );
/// Retrieve the list of edges in the adjacency matrix/graph
  void retrieveEdgeList( unsigned& nedge, std::vector<std::pair<unsigned,unsigned> >& edge_list );
///
  void getMatrixIndices( const unsigned& code, unsigned& i, unsigned& j ) const ;
/// Can we think of the matrix as an undirected graph
  bool undirectedGraph() const ;
/// Is the matrix symmetric
  bool isSymmetric() const ;
/// Get the number of rows
  unsigned getNumberOfRows() const ;
/// Get the number of columns
  unsigned getNumberOfColumns() const ;
/// Are these two nodes connected
  bool nodesAreConnected( const unsigned& iatom, const unsigned& jatom ) const ;
/// Get the cutoff that we are using for connections
  double getCutoffForConnection() const ;
///
  Vector getNodePosition( const unsigned& taskIndex ) const ;
};

}
}
#endif

