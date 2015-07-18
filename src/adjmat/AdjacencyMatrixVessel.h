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
/// This is used to keep track of what is calculated where
  std::vector<int> colvar_label;
/// The multicolvars from which we construct these quantities
  std::vector<multicolvar::MultiColvarBase*> mybasemulticolvars;
/// The vessels in these multicolvars in which the data is stored
  std::vector<vesselbase::StoreDataVessel*> mybasedata;
/// Pointer to underlying action
  AdjacencyMatrixBase* function;
/// Has the vessel been finished
  bool finished;
/// This ensures that the data is stored by the underlying multicolvars
  void buildRemoteDataStashes( const double& wtol );
/// Converts an index to the local index for that multicolvar
  unsigned convertToLocalIndex( const unsigned& index, const unsigned& mcv_code ) const ;
/// Is the underlying multicolvar active
  bool isCurrentlyActive( const unsigned& code );
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit AdjacencyMatrixVessel( const vesselbase::VesselOptions& );
/// Get the underlying adjacency matrix action object
  AdjacencyMatrixBase* getMatrixAction();
/// Ensures we use less memory for buffer in final loop
  void setBufferStart( unsigned& start );
/// Ensures that finish is set properly
  void prepare();
/// Get the position of an atom for a multicolvar
  Vector getCentralAtomPos( const unsigned& iatom ) const ;
/// Set the finished flag true
  void setFinishedTrue();
/// An overwrite of calculate to stop this being done more than once
  bool calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_index ) const ;
/// Finish the calculation
  void finish( const std::vector<double>& buffer );
/// Get the adjacency matrix
  void retrieveMatrix( DynamicList<unsigned>& myactive_elements, Matrix<double>& mymatrix );
/// Get the neighbour list based on the adjacency matrix
  void retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list );
///
  void getMatrixIndices( const unsigned& code, unsigned& i, unsigned& j );
/// Get the number of base colvars
  unsigned getNumberOfBaseColvars() const ;
///
  multicolvar::MultiColvarBase* getBaseMultiColvar( const unsigned& icolv );
///
  vesselbase::StoreDataVessel* getBaseData( const unsigned& icolv );
};

inline
void AdjacencyMatrixVessel::setBufferStart( unsigned& start ){
  if( finished ){ bufstart=start; }
  else { Vessel::setBufferStart( start ); } 
}

inline
unsigned AdjacencyMatrixVessel::convertToLocalIndex( const unsigned& index, const unsigned& mcv_code ) const {
  unsigned t1 = index;
  for(unsigned k=0;k<mcv_code;++k) t1 -= mybasemulticolvars[k]->getFullNumberOfTasks();
  return t1;
}

inline
Vector AdjacencyMatrixVessel::getCentralAtomPos( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<colvar_label.size() && colvar_label[iatom]>=0 ); 
  unsigned mmc=colvar_label[ iatom ];
  return mybasemulticolvars[mmc]->getCentralAtomPos( convertToLocalIndex(iatom,mmc) );
}

inline
bool AdjacencyMatrixVessel::isCurrentlyActive( const unsigned& code ){
  plumed_dbg_assert( code<colvar_label.size() ); unsigned mmc=colvar_label[code];
  if( mmc<0 ) return true;
  return mybasedata[mmc]->storedValueIsActive( convertToLocalIndex(code,mmc) );
}

inline
unsigned AdjacencyMatrixVessel::getNumberOfBaseColvars() const {
  return mybasemulticolvars.size();
}

inline
multicolvar::MultiColvarBase* AdjacencyMatrixVessel::getBaseMultiColvar( const unsigned& icolv ){
  return mybasemulticolvars[icolv];
}

inline
vesselbase::StoreDataVessel* AdjacencyMatrixVessel::getBaseData( const unsigned& icolv ){
  return mybasedata[icolv];
}

}
}
#endif

