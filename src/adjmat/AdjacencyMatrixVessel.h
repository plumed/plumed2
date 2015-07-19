/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#ifndef __PLUMED_multicolvar_AdjacencyMatrixVessel_h
#define __PLUMED_multicolvar_AdjacencyMatrixVessel_h

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
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit AdjacencyMatrixVessel( const vesselbase::VesselOptions& );
/// Get the underlying adjacency matrix action object
  AdjacencyMatrixBase* getMatrixAction();
/// Get the number of elements in the matrix
  unsigned getNumberOfStoredValues() const ;
/// Get the index we are storing this data inside
  unsigned getStoreIndex( const unsigned& ) const ;
/// Finish the calculation
  void finish( const std::vector<double>& buffer );
/// Get the adjacency matrix
  void retrieveMatrix( DynamicList<unsigned>& myactive_elements, Matrix<double>& mymatrix );
/// Get the neighbour list based on the adjacency matrix
  void retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list );
///
  void getMatrixIndices( const unsigned& code, unsigned& i, unsigned& j ) const ;
/// Recalculate a stored quantity
  void recalculateStoredQuantity( const unsigned& myelm, MultiValue& myvals );
};

}
}
#endif

