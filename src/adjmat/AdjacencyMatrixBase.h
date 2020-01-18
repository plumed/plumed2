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
#ifndef __PLUMED_adjmat_AdjacencyMatrixBase_h
#define __PLUMED_adjmat_AdjacencyMatrixBase_h

#include "multicolvar/MultiColvarBase.h"
#include "AdjacencyMatrixVessel.h"

namespace PLMD {
namespace adjmat {

class AdjacencyMatrixBase : public multicolvar::MultiColvarBase {
  friend class AdjacencyMatrixVessel;
  friend class ActionWithInputMatrix;
  friend class MatrixColumnSums;
  friend class MatrixRowSums;
private:
/// Used for read in of multiple connection descriptors
  unsigned connect_id;
/// Do we need to separate out the tasks for the third atoms
  bool no_third_dim_accum;
/// This is the vessel that stores the adjacency matrix
  AdjacencyMatrixVessel* mat;
/// This is used within AdjacencyMatrixVessel to recalculate matrix elements
/// whcih is useful when we are operating with lowmem
  void recalculateMatrixElement( const unsigned& myelem, MultiValue& myvals );
/// Finish the setup of the matrix
  void finishMatrixSetup( const bool& symmetric, const std::vector<AtomNumber>& all_atoms );
protected:
/// Read in a matrix involving a maximum of two species
  void readMaxTwoSpeciesMatrix( const std::string& key0, const std::string& key1, const std::string& key2, const bool& symmetric );
/// Read in a matrix involving a maximum of three species
  void readMaxThreeSpeciesMatrix( const std::string& key0, const std::string& key1, const std::string& key2, const std::string& keym, const bool& symmetric );
/// Get the dimensions of the matrix of types
  void retrieveTypeDimensions( unsigned& nrows, unsigned& ncols, unsigned& ntype ) const ;
/// Retrieve the vessel that holds the adjacency matrix
  AdjacencyMatrixVessel* getAdjacencyVessel();
/// Put the indices of the matrix elements in current atoms
  void setMatrixIndexesForTask( const unsigned& ii );
/// Add derivatives to a matrix element
  void addDerivativesOnMatrixElement( const unsigned& ielem, const unsigned& jrow, const double& df, Matrix<double>& der );
/// Read in the information on the connectors
  void parseConnectionDescriptions( const std::string& key, const bool& multiple, const unsigned& nrow_t );
protected:
/// Get the number of nodes of different types
  unsigned getNumberOfNodeTypes() const ;
/// Get the size of the vectors that were stored in the base colvars
  unsigned getSizeOfInputVectors() const ;
/// Return the group this atom is a part of
  unsigned getBaseColvarNumber( const unsigned& ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit AdjacencyMatrixBase(const ActionOptions&);
/// Create the connection object
  virtual void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) = 0;
/// None of these things are allowed
  bool isPeriodic() { return false; }
  Vector getCentralAtom() { plumed_merror("cannot find central atoms for adjacency matrix actions"); Vector dum; return dum; }
/// Get the atom number
  AtomNumber getAbsoluteIndexOfCentralAtom( const unsigned& i ) const ;
};

inline
AdjacencyMatrixVessel* AdjacencyMatrixBase::getAdjacencyVessel() {
  return mat;
}

inline
unsigned AdjacencyMatrixBase::getBaseColvarNumber( const unsigned& inum ) const {
  if( atom_lab[inum].first>0 ) return atom_lab[inum].first-1;
  return 0;
}

inline
AtomNumber AdjacencyMatrixBase::getAbsoluteIndexOfCentralAtom( const unsigned& iatom ) const {
  if( atom_lab[iatom].first>0 ) {
    unsigned mmc=atom_lab[ iatom ].first - 1;
    return mybasemulticolvars[mmc]->getAbsoluteIndexOfCentralAtom( atom_lab[iatom].second );
  }
  return ActionAtomistic::getAbsoluteIndex( atom_lab[iatom].second );
}

}
}

#endif
