/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#ifndef __PLUMED_tools_MultiValue_h
#define __PLUMED_tools_MultiValue_h

#include "Exception.h"
#include "Vector.h"
#include "Tensor.h"
#include <vector>
#include <cstddef>

namespace PLMD {

class MultiValue {
  friend class ActionWithVector;
private:
/// The index of the task we are currently performing
  std::size_t task_index, task2_index;
/// Values of quantities
  std::vector<double> values;
/// Number of derivatives per value
  unsigned nderivatives;
/// Derivatives
  std::vector<double> derivatives;
/// Matrix asserting which values have derivatives
  std::vector<bool> hasderiv;
/// Tempory value
  double tmpval;
/// Lists of active variables
  std::vector<unsigned> nactive, active_list;
/// Tempory vector of derivatives (used for calculating quotients
  std::vector<double> tmpder;
/// Logical to check if any derivatives were set
  bool atLeastOneSet;
/// Are we in this for a call on vectors
  bool vector_call;
  unsigned nindices, nsplit;
/// This allows us to store matrix elements
  unsigned nmatrix_cols;
  std::vector<double> matrix_row_stash;
  std::vector<double> matrix_force_stash;
  std::vector<unsigned> matrix_bookeeping;
/// These are used to store the indices that have derivatives wrt to at least one
/// of the elements in a matrix
  std::vector<unsigned> matrix_row_nderivatives;
  std::vector<std::vector<unsigned> > matrix_row_derivative_indices;
/// This is a fudge to save on vector resizing in MultiColvar
  std::vector<unsigned> indices;
  std::vector<Vector> tmp_atoms;
  std::vector<std::vector<Vector> > tmp_atom_der;
  std::vector<Tensor> tmp_atom_virial;
  std::vector<std::vector<double> > tmp_vectors;
public:
  MultiValue( const std::size_t& nvals, const std::size_t& nder, const std::size_t& nmat=0, const std::size_t& maxcol=0, const std::size_t& nbook=0 );
  void resize( const std::size_t& nvals, const std::size_t& nder, const std::size_t& nmat=0, const std::size_t& maxcol=0, const std::size_t& nbook=0 );
/// Set the task index prior to the loop
  void setTaskIndex( const std::size_t& tindex );
///
  std::size_t getTaskIndex() const ;
///
  void setSecondTaskIndex( const std::size_t& tindex );
/// Get the task index
  std::size_t getSecondTaskIndex() const ;
///
  void setSplitIndex( const std::size_t& nat );
  std::size_t getSplitIndex() const ;
///
  void setNumberOfIndices( const std::size_t& nat );
  std::size_t getNumberOfIndices() const ;
///
  std::vector<unsigned>& getIndices();
  std::vector<Vector>& getAtomVector();
/// Get the number of values in the stash
  unsigned getNumberOfValues() const ;
/// Get the number of derivatives in the stash
  unsigned getNumberOfDerivatives() const ;
/// Get references to some memory. These vectors allow us to
/// avoid doing lots of resizing of vectors in MultiColvarTemplate
  std::vector<Vector>& getFirstAtomVector();
  std::vector<std::vector<Vector> >& getFirstAtomDerivativeVector();
  const std::vector<std::vector<Vector> >& getConstFirstAtomDerivativeVector() const ;
  std::vector<Tensor>& getFirstAtomVirialVector();
  void resizeTemporyVector(const unsigned& n );
  std::vector<double>& getTemporyVector(const unsigned& ind );
///
  bool inVectorCall() const ;
/// Set value numbered
  void setValue( const std::size_t&,  const double& );
/// Add value numbered
  void addValue( const std::size_t&,  const double& );
/// Add derivative
  void addDerivative( const std::size_t&, const std::size_t&, const double& );
/// Set the value of the derivative
  void setDerivative( const std::size_t& ival, const std::size_t& jder, const double& der);
/// Return the ith value
  double get( const std::size_t& ) const ;
/// Return a derivative value
  double getDerivative( const std::size_t&, const std::size_t& ) const ;
/// Clear all values
  void clearAll();
/// Clear the derivatives
  void clearDerivatives( const unsigned& );
/// Clear a value
  void clear( const unsigned& );
/// Functions for accessing active list
  bool updateComplete();
  void emptyActiveMembers();
  void putIndexInActiveArray( const unsigned & );
  void updateIndex( const unsigned& );
///
  void updateIndex( const std::size_t&, const std::size_t& );
///
  unsigned getActiveIndex( const std::size_t&, const std::size_t& ) const ;
///
  void clearActiveMembers( const std::size_t& ival );
///
  unsigned getNumberActive( const std::size_t& ival ) const ;
///
  unsigned getActiveIndex( const unsigned& ) const ;
/// Get the matrix bookeeping array
  const std::vector<unsigned> & getMatrixBookeeping() const ;
  void stashMatrixElement( const unsigned& nmat, const unsigned& rowstart, const unsigned& jcol, const double& val );
  double getStashedMatrixElement( const unsigned& nmat, const unsigned& jcol ) const ;
/// Get the bookeeping stuff for the derivatives wrt to rows of matrix
  void setNumberOfMatrixRowDerivatives( const unsigned& nmat, const unsigned& nind );
  unsigned getNumberOfMatrixRowDerivatives( const unsigned& nmat ) const ;
  std::vector<unsigned>& getMatrixRowDerivativeIndices( const unsigned& nmat );
/// Stash the forces on the matrix
  void addMatrixForce( const unsigned& imat, const unsigned& jind, const double& f );
  double getStashedMatrixForce( const unsigned& imat, const unsigned& jind ) const ;
};

inline
unsigned MultiValue::getNumberOfValues() const {
  return values.size();
}

inline
unsigned MultiValue::getNumberOfDerivatives() const {
  return nderivatives; //derivatives.ncols();
}

inline
double MultiValue::get( const std::size_t& ival ) const {
  plumed_dbg_assert( ival<=values.size() );
  return values[ival];
}

inline
void MultiValue::setValue( const std::size_t& ival,  const double& val) {
  plumed_dbg_assert( ival<=values.size() );
  values[ival]=val;
}

inline
void MultiValue::addValue( const std::size_t& ival,  const double& val) {
  plumed_dbg_assert( ival<=values.size() );
  values[ival]+=val;
}

inline
void MultiValue::addDerivative( const std::size_t& ival, const std::size_t& jder, const double& der) {
  plumed_dbg_assert( ival<=values.size() && jder<nderivatives ); atLeastOneSet=true;
  hasderiv[nderivatives*ival+jder]=true; derivatives[nderivatives*ival+jder] += der;
}

inline
void MultiValue::setDerivative( const std::size_t& ival, const std::size_t& jder, const double& der) {
  plumed_dbg_assert( ival<=values.size() && jder<nderivatives ); atLeastOneSet=true;
  hasderiv[nderivatives*ival+jder]=true; derivatives[nderivatives*ival+jder]=der;
}


inline
double MultiValue::getDerivative( const std::size_t& ival, const std::size_t& jder ) const {
  plumed_dbg_assert( ival<values.size() && jder<nderivatives );
  return derivatives[nderivatives*ival+jder];
}

inline
void MultiValue::updateIndex( const std::size_t& ival, const std::size_t& jder ) {
  plumed_dbg_assert( ival<values.size() && jder<nderivatives );
#ifdef DNDEBUG
  for(unsigned i=0; i<nactive[ival]; ++i) plumed_dbg_assert( active_list[nderivatives*ival+nactive[ival]]!=jder );
#endif
  if( hasderiv[nderivatives*ival+jder] ) {
    plumed_dbg_assert( nactive[ival]<nderivatives);
    active_list[nderivatives*ival+nactive[ival]]=jder;
    nactive[ival]++;
  }
}

inline
unsigned MultiValue::getNumberActive( const std::size_t& ival ) const {
  plumed_dbg_assert( ival<nactive.size() );
  return nactive[ival];
}

inline
unsigned MultiValue::getActiveIndex( const std::size_t& ival, const std::size_t& ind ) const {
  plumed_dbg_assert( ind<nactive[ival] );
  return active_list[nderivatives*ival+ind];
}

inline
void MultiValue::setSplitIndex( const std::size_t& nat ) {
  nsplit = nat;
}

inline
std::size_t MultiValue::getSplitIndex() const {
  return nsplit;
}

inline
void MultiValue::setNumberOfIndices( const std::size_t& nat ) {
  nindices = nat;
}

inline
std::size_t MultiValue::getNumberOfIndices() const {
  return nindices;
}


inline
bool MultiValue::inVectorCall() const {
  return (matrix_row_nderivatives.size()>0 && vector_call);
}

inline
void MultiValue::clearActiveMembers( const std::size_t& ival ) {
  nactive[ival]=0;
}

inline
void MultiValue::setTaskIndex( const std::size_t& tindex ) {
  task_index = tindex;
}

inline
std::size_t MultiValue::getTaskIndex() const {
  return task_index;
}

inline
void MultiValue::setSecondTaskIndex( const std::size_t& tindex ) {
  task2_index = tindex;
}

inline
std::size_t MultiValue::getSecondTaskIndex() const {
  return task2_index;
}

inline
std::vector<unsigned>& MultiValue::getIndices() {
  return indices;
}

inline
std::vector<Vector>& MultiValue::getAtomVector() {
  return tmp_atoms;
}

inline
std::vector<Vector>& MultiValue::getFirstAtomVector() {
  return tmp_atoms;
}

inline
std::vector<std::vector<Vector> >& MultiValue::getFirstAtomDerivativeVector() {
  return tmp_atom_der;
}

inline
const std::vector<std::vector<Vector> >& MultiValue::getConstFirstAtomDerivativeVector() const {
  return tmp_atom_der;
}

inline
std::vector<Tensor>& MultiValue::getFirstAtomVirialVector() {
  return tmp_atom_virial;
}

inline
void MultiValue::stashMatrixElement( const unsigned& nmat, const unsigned& rowstart, const unsigned& jcol, const double& val ) {
  plumed_dbg_assert( jcol<nmatrix_cols && rowstart + matrix_bookeeping[rowstart]<matrix_bookeeping.size() && nmatrix_cols*nmat + matrix_bookeeping[rowstart]<matrix_row_stash.size() );
  matrix_bookeeping[rowstart]++; matrix_bookeeping[rowstart + matrix_bookeeping[rowstart]]=jcol; matrix_row_stash[ nmatrix_cols*nmat + jcol] = val;
}

inline
double MultiValue::getStashedMatrixElement( const unsigned& nmat, const unsigned& jcol ) const {
  plumed_dbg_assert( nmatrix_cols*nmat + jcol<matrix_row_stash.size() );
  return matrix_row_stash[ nmatrix_cols*nmat + jcol ];
}

inline
const std::vector<unsigned> & MultiValue::getMatrixBookeeping() const {
  return matrix_bookeeping;
}

inline
void MultiValue::setNumberOfMatrixRowDerivatives( const unsigned& nmat, const unsigned& nind ) {
  plumed_dbg_assert( nmat<matrix_row_nderivatives.size() && nind<=matrix_row_derivative_indices[nmat].size() );
  matrix_row_nderivatives[nmat]=nind;
}

inline
unsigned MultiValue::getNumberOfMatrixRowDerivatives( const unsigned& nmat ) const {
  plumed_dbg_assert( nmat<matrix_row_nderivatives.size() ); return matrix_row_nderivatives[nmat];
}

inline
std::vector<unsigned>& MultiValue::getMatrixRowDerivativeIndices( const unsigned& nmat ) {
  plumed_dbg_assert( nmat<matrix_row_nderivatives.size() ); return matrix_row_derivative_indices[nmat];
}

inline
void MultiValue::addMatrixForce( const unsigned& imat, const unsigned& jind, const double& f ) {
  matrix_force_stash[imat*nderivatives + jind]+=f;
}

inline
double MultiValue::getStashedMatrixForce( const unsigned& imat, const unsigned& jind ) const {
  return matrix_force_stash[imat*nderivatives + jind];
}

inline
void MultiValue::resizeTemporyVector(const unsigned& n ) {
  if( n>tmp_vectors.size() ) tmp_vectors.resize(n);
}

inline
std::vector<double>& MultiValue::getTemporyVector(const unsigned& ind ) {
  plumed_dbg_assert( ind<tmp_vectors.size() );
  return tmp_vectors[ind];
}

}
#endif
