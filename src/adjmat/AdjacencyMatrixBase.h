/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
friend class MatrixSummationBase;
private:
/// Used for read in of multiple connection descriptors
  unsigned connect_id;
/// The tolerance to use to decide whether or not to incorporate atoms
  double wtolerance;
/// Do we need to separate out the tasks for the third atoms
  bool no_third_dim_accum;
/// This is the vessel that stores the adjacency matrix
  AdjacencyMatrixVessel* mat;
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
  bool parseAtomList(const std::string& key, const int& num, std::vector<AtomNumber>& t);
/// Get the number of nodes of different types
  unsigned getNumberOfNodeTypes() const ;
/// Get the size of the vectors that were stored in the base colvars
  unsigned getSizeOfInputVectors() const ;
/// Request the atoms
  void requestAtoms( const std::vector<AtomNumber>& atoms, const bool& symmetric, const bool& true_square, const std::vector<unsigned>& dims );
/// Return the group this atom is a part of
  unsigned getBaseColvarNumber( const unsigned& ) const ;
/// Add some derivatives to the relevant orientation
  void addOrientationDerivatives( const unsigned&, const unsigned& , const std::vector<double>& , multicolvar::AtomValuePack& ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit AdjacencyMatrixBase(const ActionOptions&);
/// Create the connection object
  virtual void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::string& desc ) = 0;
/// None of these things are allowed
  bool isPeriodic(){ return false; }
  Vector getCentralAtom(){ plumed_merror("cannot find central atoms for adjacency matrix actions"); Vector dum; return dum; }
/// Get the absolute index of an atom
//  AtomNumber getAbsoluteIndexOfCentralAtom(const unsigned& i) const ;
/// Transforms the stored values in whatever way is required
  virtual double transformStoredValues( const std::vector<double>& myvals, unsigned& vout, double& df ) const ;
/// Used to check for connections between atoms
  virtual bool checkForConnection( const std::vector<double>& myvals ) const=0;
/// Get the atom number
  AtomNumber getAbsoluteIndexOfCentralAtom( const unsigned& i ) const ; 
};

inline
AdjacencyMatrixVessel* AdjacencyMatrixBase::getAdjacencyVessel(){
  return mat;
}

inline
unsigned AdjacencyMatrixBase::getBaseColvarNumber( const unsigned& inum ) const {
  if( inum<colvar_label.size() ) return colvar_label[inum]; 
  return 0;
}

inline
AtomNumber AdjacencyMatrixBase::getAbsoluteIndexOfCentralAtom( const unsigned& iatom ) const {
  if( iatom<colvar_label.size() ){
      unsigned mmc=colvar_label[ iatom ];
      return mybasemulticolvars[mmc]->getAbsoluteIndexOfCentralAtom( convertToLocalIndex(iatom,mmc) );
  }
  return ActionAtomistic::getAbsoluteIndex( iatom );
}


inline 
void AdjacencyMatrixBase::getOrientationVector( const unsigned& ind, const bool& normed, std::vector<double>& orient ) const {
  plumed_dbg_assert( ind<colvar_label.size() ); unsigned mmc=colvar_label[ind];
  plumed_assert( !mybasemulticolvars[mmc]->weightWithDerivatives() ); plumed_dbg_assert( mybasedata[mmc]->storedValueIsActive( convertToLocalIndex(ind,mmc) ) );
  mybasedata[mmc]->retrieveValueWithIndex( convertToLocalIndex(ind,mmc), normed, orient );
}

inline
double AdjacencyMatrixBase::transformStoredValues( const std::vector<double>& myvals, unsigned& vout, double& df  ) const {
  plumed_dbg_assert( myvals.size()==2 ); vout=1; df=1; return myvals[1]; 
}

}
}

#endif
