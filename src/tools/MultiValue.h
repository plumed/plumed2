/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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

#include <vector>
#include "Exception.h"
#include "Vector.h"

namespace PLMD {

class MultiValue {
friend class ActionWithValue;
private:
/// The index of the task we are currently performing
  unsigned task_index;
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
  unsigned nmatrix_cols;
  std::vector<double> matrix_element_stash;
/// This is a fudge to save on vector resizing in MultiColvar
  bool vector_call;
  unsigned nindices;
  std::vector<unsigned> indices, sort_indices;
  std::vector<Vector> tmp_atoms;
public:
  MultiValue( const unsigned& nvals, const unsigned& nder, const unsigned ncols=0, const unsigned nmat=0 );
  void resize( const unsigned&, const unsigned& );
/// Set the task index prior to the loop
  void setTaskIndex( const unsigned& tindex );
/// Get the task index
  unsigned getTaskIndex() const ;
///
  void setNumberOfIndices( const unsigned& nat );
  unsigned getNumberOfIndices() const ;
  std::vector<unsigned>& getIndices();
  std::vector<unsigned>& getSortIndices();
  std::vector<Vector>& getAtomVector();
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
/// Clear all values
  void clearAll();
/// Clear a value
  void clear( const unsigned& );
  void updateIndex( const unsigned& , const unsigned& );
///
  unsigned getNumberActive( const unsigned& ) const ;
///
  unsigned getActiveIndex( const unsigned& , const unsigned& ) const ;
/// 
  void stashMatrixElement( const unsigned& imat, const unsigned& jind, const double& val );
///
  double getStashedMatrixElement( const unsigned& imat, const unsigned& jind ) const ;
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
  plumed_dbg_assert( ival<values.size() && jder<nderivatives && hasderiv[nderivatives*ival+jder] );
  return derivatives[nderivatives*ival+jder];
}

inline
void MultiValue::updateIndex( const unsigned& ival, const unsigned& jder ) {
  plumed_dbg_assert( ival<values.size() && jder<nderivatives );
#ifdef DNDEBUG
  for(unsigned i=0;i<nactive[ival];++i) plumed_dbg_assert( active_list[nderivatives*ival+nactive[ival]]!=jder ); 
#endif
  if( hasderiv[nderivatives*ival+jder] ){ active_list[nderivatives*ival+nactive[ival]]=jder; nactive[ival]++; } 
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
std::vector<unsigned>& MultiValue::getSortIndices() {
  return sort_indices;
}

inline
std::vector<Vector>& MultiValue::getAtomVector() {
  return tmp_atoms;
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
void MultiValue::stashMatrixElement( const unsigned& imat, const unsigned& jind, const double& val ) {
  matrix_element_stash[imat*nmatrix_cols + jind] = val;
}

inline
double MultiValue::getStashedMatrixElement( const unsigned& imat, const unsigned& jind ) const {
  return matrix_element_stash[imat*nmatrix_cols + jind];
}

inline
bool MultiValue::inVectorCall() const {
  return vector_call;
}

}
#endif
