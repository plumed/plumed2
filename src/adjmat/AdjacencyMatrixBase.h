/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#ifndef __PLUMED_adjmat_AdjacencyMatrixBase_h
#define __PLUMED_adjmat_AdjacencyMatrixBase_h

#include "multicolvar/MultiColvarBase.h"
#include "multicolvar/InputMultiColvarSet.h"
#include "AdjacencyMatrixVessel.h"

namespace PLMD {
namespace adjmat {

class AdjacencyMatrixBase : public multicolvar::MultiColvarBase {
friend class AdjacencyMatrixVessel;
friend class ActionWithInputMatrix;
friend class MatrixSummationBase;
private:
/// Used for read in of multiple connection descriptors
  unsigned connect_id;
/// The tolerance to use to decide whether or not to incorporate atoms
  double wtolerance;
/// This is the vessel that stores the adjacency matrix
  AdjacencyMatrixVessel* mat;
/// This stores the base colvars
  multicolvar::InputMultiColvarSet myinputdata;
/// This is used within AdjacencyMatrixVessel to recalculate matrix elements
/// whcih is useful when we are operating with lowmem
  void recalculateMatrixElement( const unsigned& myelem, MultiValue& myvals );
protected:
/// Retrieve the vessel that holds the adjacency matrix
  AdjacencyMatrixVessel* getAdjacencyVessel();
/// Get the vector for a particular node
  void getOrientationVector( const unsigned& ind, const bool& normed, std::vector<double>& orient0 ) const ; 
/// Put the indices of the matrix elements in current atoms
  void setMatrixIndexesForTask( const unsigned& ii );
/// Add derivatives to a matrix element
  void addDerivativesOnMatrixElement( const unsigned& ielem, const unsigned& jrow, const double& df, Matrix<double>& der );
/// Read in the information on the connectors
  void parseConnectionDescriptions( const std::string& key, const unsigned& nrow_t );
protected:
/// Read the list of atoms involved in this colvar
  void parseAtomList(const std::string& key, const int& num, const bool& isnodes, std::vector<AtomNumber>& t);
/// Get the number of nodes of different types
  unsigned getNumberOfNodeTypes() const ;
/// Get the total number of nodes
  unsigned getNumberOfNodes() const ;
/// Get the size of the vectors that were stored in the base colvars
  unsigned getSizeOfInputVectors() const ;
/// Request the atoms
  void requestAtoms( const std::vector<AtomNumber>& atoms, const bool& symmetric, const unsigned& nrows );
/// Return the group this atom is a part of
  unsigned getBaseColvarNumber( const unsigned& ) const ;
/// Add some derivatives to the relevant atom
  void addAtomDerivatives( const unsigned& , const Vector& , multicolvar::AtomValuePack& ) const ;
/// Add some derivatives to the relevant orientation
  void addOrientationDerivatives( const unsigned& iatom, const std::vector<double>& der, multicolvar::AtomValuePack& myatoms ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit AdjacencyMatrixBase(const ActionOptions&);
/// This returns the position of an atom for link cells 
  Vector getPositionOfAtomForLinkCells( const unsigned& iatom ) const ;
/// Tells us if a particular index is active 
  bool isCurrentlyActive( const unsigned& bno, const unsigned& code );
/// Update the list of active atoms in the underlying atom value pack 
  void updateActiveAtoms( multicolvar::AtomValuePack& myatoms ) const ;
/// Calculation routine
  void calculate();
/// Create the connection object
  virtual void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::string& desc ) = 0;
/// None of these things are allowed
  bool isPeriodic(){ return false; }
  Vector getCentralAtom(){ plumed_merror("cannot find central atoms for adjacency matrix actions"); Vector dum; return dum; }
};

inline
AdjacencyMatrixVessel* AdjacencyMatrixBase::getAdjacencyVessel(){
  return mat;
}

inline
unsigned AdjacencyMatrixBase::getBaseColvarNumber( const unsigned& inum ) const {
  return myinputdata.getBaseColvarNumber(inum);
}

inline 
void AdjacencyMatrixBase::getOrientationVector( const unsigned& ind, const bool& normed, std::vector<double>& orient ) const {
  myinputdata.getVectorForTask( ind, normed, orient );
}

}
}

#endif
