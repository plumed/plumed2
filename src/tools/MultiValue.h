/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include <vector>
#include <cstddef>

namespace PLMD {

class MultiValue {
  friend class ActionWithValue;
private:
/// The index of the task we are currently performing
  unsigned task_index, task2_index;
/// Values of quantities
  std::vector<double> values;
/// Number of derivatives per value
  unsigned nderivatives;
/// Derivatives
  std::vector<double> derivatives;
/// Matrix asserting which values have derivatives
  std::vector<bool> hasderiv;
/// Lists of active variables
  std::vector<unsigned> nactive, active_list;
/// Logical to check if any derivatives were set
  bool atLeastOneSet;
/// This allows us to store matrix elements
  unsigned nmatrix_cols, nmat_force;
  bool rerunning_matrix;
  std::vector<unsigned> matrix_element_nind;
  std::vector<unsigned> matrix_element_indices;
  std::vector<double> matrix_element_stash;
  std::vector<double> matrix_force_stash;
/// This is a fudge to save on vector resizing in MultiColvar
  bool vector_call;
  unsigned nindices, nfblock, nsplit;
  std::vector<unsigned> indices;
  std::vector<unsigned> mat_nindices;
  std::vector<std::vector<unsigned> > mat_indices;
  std::vector<std::vector<Vector> > tmp_atoms;
/// Tempory storage space for doing stuff with symmetry functions
  unsigned symfunc_index;
  std::vector<std::vector<double> > symfunc_tmp_derivs;
public:
  MultiValue( const std::size_t& nvals, const std::size_t& nder, const std::size_t ncols=0, const std::size_t nmat=0, const std::size_t nfder=0 );
  void resize( const std::size_t&, const std::size_t&, const std::size_t&, const std::size_t& );
/// Set the task index prior to the loop
  void setTaskIndex( const unsigned& tindex );
/// Get the task index
  unsigned getTaskIndex() const ;
///
  void setSecondTaskIndex( const unsigned& tindex );
/// Get the task index
  unsigned getSecondTaskIndex() const ;
/// Tempory indices
  void setSplitIndex( const unsigned& nat );
  unsigned getSplitIndex() const ;
  void setNumberOfIndicesInFirstBlock( const unsigned& nat );
  unsigned getNumberOfIndicesInFirstBlock( ) const ;
  void setNumberOfIndices( const unsigned& nat );
  unsigned getNumberOfIndices() const ;
  std::vector<unsigned>& getIndices();
  const std::vector<unsigned>& getIndices() const ;
/// Tempory matrix indices
  void setNumberOfMatrixIndices( const unsigned& nmat, const unsigned& nind );
  unsigned getNumberOfMatrices() const ;
  unsigned getNumberOfColumns() const ;
  unsigned getNumberOfMatrixIndices( const unsigned& nmat ) const ;
  std::vector<unsigned>& getMatrixIndices( const unsigned& nmat );
  std::vector<Vector>& getFirstAtomVector();
  const std::vector<Vector>& getAtomVector() const ;
  std::vector<Vector>& getSecondAtomVector();
/// Get the number of values in the stash
  unsigned getNumberOfValues() const ;
/// Get the number of derivatives in the stash
  unsigned getNumberOfDerivatives() const ;
///
  bool inVectorCall() const ;
/// Set value numbered
  void setValue( const unsigned&,  const double& );
/// Add value numbered
  void addValue( const unsigned&,  const double& );
/// Add derivative
  void addDerivative( const unsigned&, const unsigned&, const double& );
/// Set the value of the derivative
  void setDerivative( const unsigned& ival, const unsigned& jder, const double& der);
/// Return the ith value
  double get( const unsigned& ) const ;
/// Return a derivative value
  double getDerivative( const unsigned&, const unsigned& ) const ;
/// Clear bookeeping arrays for matrix stuff
  void clearMatrixBookeepingArrays();
/// Clear all values
  void clearAll();
/// Clear a value
  void clear( const unsigned& );
///
  void clearActiveMembers( const unsigned& vv );
  void updateIndex( const unsigned&, const unsigned& );
///
  unsigned getNumberActive( const unsigned& ) const ;
///
  unsigned getActiveIndex( const unsigned&, const unsigned& ) const ;
///
  void stashMatrixElement( const unsigned& imat, const unsigned& jind, const double& val );
///
  unsigned getNumberOfStashedMatrixElements( const unsigned& imat ) const ;
///
  unsigned getStashedMatrixIndex( const unsigned& imat, const unsigned& jind ) const ;
///
  double getStashedMatrixElement( const unsigned& imat, const unsigned& jind ) const ;
///
  void addMatrixForce( const unsigned& imat, const unsigned& jind, const double& f );
///
  double getStashedMatrixForce( const unsigned& imat, const unsigned& jind ) const ;
///
  void clearStoredForces();
///
  void setMatrixStashForRerun();
///
  bool inMatrixRerun() const ;
///
  void setMatrixStashForNormalRun();
///
  void setSymfuncTemporyIndex( const unsigned& ind );
///
  unsigned getSymfuncTemporyIndex() const ;
///
  std::vector<double>& getSymfuncTemporyDerivatives( const unsigned& ind );
};

inline
unsigned MultiValue::getNumberOfValues() const {
  return values.size();
}

inline
unsigned MultiValue::getNumberOfDerivatives() const {
  return nderivatives;
}

inline
double MultiValue::get( const unsigned& ival ) const {
  plumed_dbg_assert( ival<values.size() );
  return values[ival];
}

inline
void MultiValue::setValue( const unsigned& ival,  const double& val) {
  plumed_dbg_assert( ival<values.size() );
  values[ival]=val;
}

inline
void MultiValue::addValue( const unsigned& ival,  const double& val) {
  plumed_dbg_assert( ival<values.size() );
  values[ival]+=val;
}

inline
void MultiValue::addDerivative( const unsigned& ival, const unsigned& jder, const double& der) {
  plumed_dbg_assert( ival<values.size() && jder<nderivatives ); atLeastOneSet=true;
  hasderiv[nderivatives*ival+jder]=true; derivatives[nderivatives*ival+jder] += der;
}

inline
void MultiValue::setDerivative( const unsigned& ival, const unsigned& jder, const double& der) {
  plumed_dbg_assert( ival<values.size() && jder<nderivatives ); atLeastOneSet=true;
  hasderiv[nderivatives*ival+jder]=true; derivatives[nderivatives*ival+jder]=der;
}

inline
double MultiValue::getDerivative( const unsigned& ival, const unsigned& jder ) const {
  plumed_dbg_assert( ival<values.size() && jder<nderivatives );
  if( !hasderiv[nderivatives*ival+jder] ) return 0.0;
  return derivatives[nderivatives*ival+jder];
}

inline
void MultiValue::updateIndex( const unsigned& ival, const unsigned& jder ) {
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
unsigned MultiValue::getNumberActive( const unsigned& ival ) const {
  plumed_dbg_assert( ival<nactive.size() );
  return nactive[ival];
}

inline
unsigned MultiValue::getActiveIndex( const unsigned& ival, const unsigned& ind ) const {
  plumed_dbg_assert( ind<nactive[ival] );
  return active_list[nderivatives*ival+ind];
}

inline
void MultiValue::setSplitIndex( const unsigned& nat ) {
  nsplit = nat;
}

inline
unsigned MultiValue::getSplitIndex() const {
  return nsplit;
}

inline
void MultiValue::setNumberOfIndicesInFirstBlock( const unsigned& nat ) {
  nfblock = nat;
}

inline
unsigned MultiValue::getNumberOfIndicesInFirstBlock() const {
  return nfblock;
}

inline
void MultiValue::setNumberOfIndices( const unsigned& nat ) {
  nindices = nat;
}

inline
unsigned MultiValue::getNumberOfIndices() const {
  return nindices;
}

inline
std::vector<unsigned>& MultiValue::getIndices() {
  return indices;
}

inline
const std::vector<unsigned>& MultiValue::getIndices() const {
  return indices;
}

inline
std::vector<Vector>& MultiValue::getFirstAtomVector() {
  return tmp_atoms[0];
}

inline
const std::vector<Vector>& MultiValue::getAtomVector() const {
  return tmp_atoms[0];
}

inline
std::vector<Vector>& MultiValue::getSecondAtomVector() {
  return tmp_atoms[1];
}

inline
void MultiValue::setTaskIndex( const unsigned& tindex ) {
  task_index = tindex;
}

inline
unsigned MultiValue::getTaskIndex() const {
  return task_index;
}

inline
void MultiValue::setSecondTaskIndex( const unsigned& tindex ) {
  task2_index = tindex;
}

inline
unsigned MultiValue::getSecondTaskIndex() const {
  return task2_index;
}

inline
unsigned MultiValue::getNumberOfStashedMatrixElements( const unsigned& imat ) const {
  plumed_dbg_assert( imat<matrix_element_nind.size() );
  return matrix_element_nind[imat];
}

inline
unsigned MultiValue::getStashedMatrixIndex( const unsigned& imat, const unsigned& jind ) const {
  plumed_dbg_assert( imat<matrix_element_nind.size() && jind<matrix_element_nind[imat] );
  return matrix_element_indices[imat*nmatrix_cols + jind];
}

inline
void MultiValue::stashMatrixElement( const unsigned& imat, const unsigned& jind, const double& val ) {
  if( rerunning_matrix ) return;
  plumed_dbg_assert( imat<matrix_element_nind.size() && jind<nmatrix_cols );
  matrix_element_indices[imat*nmatrix_cols + matrix_element_nind[imat]] = jind;
  matrix_element_nind[imat]++; matrix_element_stash[imat*nmatrix_cols + jind] = val;
}

inline
void MultiValue::addMatrixForce( const unsigned& imat, const unsigned& jind, const double& f ) {
  matrix_force_stash[imat*nmat_force + jind]+=f;
}

inline
double MultiValue::getStashedMatrixForce( const unsigned& imat, const unsigned& jind ) const {
  return matrix_force_stash[imat*nmat_force + jind];
}

inline
void MultiValue::clearStoredForces() {
  std::fill(matrix_force_stash.begin(),matrix_force_stash.end(),0);
}

inline
double MultiValue::getStashedMatrixElement( const unsigned& imat, const unsigned& jind ) const {
  return matrix_element_stash[imat*nmatrix_cols + jind];
}

inline
void MultiValue::setMatrixStashForRerun() {
  rerunning_matrix=true;
}

inline
void MultiValue::setMatrixStashForNormalRun() {
  rerunning_matrix=false;
}

inline
bool MultiValue::inMatrixRerun() const {
  return rerunning_matrix;
}

inline
bool MultiValue::inVectorCall() const {
  return (nmatrix_cols>0 && vector_call);
}

inline
void MultiValue::clearActiveMembers( const unsigned& ival ) {
  nactive[ival]=0;
}

inline
void MultiValue::setNumberOfMatrixIndices( const unsigned& nmat, const unsigned& nind ) {
  plumed_dbg_assert( nmat<mat_nindices.size() && nind<=mat_indices[nmat].size() );
  mat_nindices[nmat]=nind;
}

inline
unsigned MultiValue::getNumberOfMatrixIndices( const unsigned& nmat ) const {
  plumed_dbg_assert( nmat<mat_nindices.size() ); return mat_nindices[nmat];
}

inline
std::vector<unsigned>& MultiValue::getMatrixIndices( const unsigned& nmat ) {
  plumed_dbg_assert( nmat<mat_nindices.size() ); return mat_indices[nmat];
}


inline
void MultiValue::setSymfuncTemporyIndex( const unsigned& ind ) {
  symfunc_index=ind;
}

inline
unsigned MultiValue::getSymfuncTemporyIndex() const {
  return symfunc_index;
}

inline
std::vector<double>& MultiValue::getSymfuncTemporyDerivatives( const unsigned& ind ) {
  plumed_dbg_assert( ind<symfunc_tmp_derivs.size() ); return symfunc_tmp_derivs[ind];
}

inline
unsigned MultiValue::getNumberOfMatrices() const {
  return matrix_element_nind.size();
}

inline
unsigned MultiValue::getNumberOfColumns() const {
  return nmatrix_cols;
}

}
#endif
