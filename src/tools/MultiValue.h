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
#include "DynamicList.h"
#include <vector>
#include <cstddef>

namespace PLMD {

class MultiValue {
friend class ActionWithVector;
private:
/// The index of the task we are currently performing
  std::size_t task_index, task2_index;
/// Used to ensure rapid accumulation of derivatives
  DynamicList<unsigned> hasDerivatives;
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
/// This allows us to store matrix elements
  unsigned nmatrix_cols;
/// This is a fudge to save on vector resizing in MultiColvar
  std::vector<unsigned> indices, sort_indices;
  std::vector<Vector> tmp_atoms;
  std::vector<std::vector<Vector> > tmp_atom_der;
  std::vector<Tensor> tmp_atom_virial;
  std::vector<std::vector<double> > tmp_vectors;
public:
  MultiValue( const std::size_t& nvals, const std::size_t& nder, const std::size_t& ncols=0, const std::size_t& nmat=0 );
  void resize( const std::size_t& nvals, const std::size_t& nder, const std::size_t& ncols=0, const std::size_t& nmat=0 );
/// Set the task index prior to the loop
  void setTaskIndex( const std::size_t& tindex );
///
  std::size_t getTaskIndex() const ;
  std::vector<unsigned>& getIndices();
  std::vector<unsigned>& getSortIndices();
  std::vector<Vector>& getAtomVector();
/// Get the number of values in the stash
  unsigned getNumberOfValues() const ;
/// Get the number of derivatives in the stash
  unsigned getNumberOfDerivatives() const ;
/// Get references to some memory. These vectors allow us to 
/// avoid doing lots of resizing of vectors in MultiColvarTemplate 
  std::vector<Vector>& getFirstAtomVector();
  std::vector<std::vector<Vector> >& getFirstAtomDerivativeVector();
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
/// Add to the tempory value
  void addTemporyValue( const double& val );
/// Add tempory derivatives - this is used for calculating quotients
  void addTemporyDerivative( const unsigned& jder, const double& der );
/// Set the value of the derivative
  void setDerivative( const std::size_t& ival, const std::size_t& jder, const double& der);
/// Return the ith value
  double get( const std::size_t& ) const ;
/// Return a derivative value
  double getDerivative( const std::size_t&, const std::size_t& ) const ;
/// Get one of the tempory derivatives
  double getTemporyDerivative( const unsigned& jder ) const ;
/// Clear all values
  void clearAll();
/// Clear the tempory derivatives
  void clearTemporyDerivatives();
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
  unsigned getNumberActive( const std::size_t& ) const ;
///
  unsigned getActiveIndex( const std::size_t&, const std::size_t& ) const ;  
///
  void clearActiveMembers( const std::size_t& ival );
  void sortActiveList();
  void completeUpdate();
  void updateDynamicList();
  bool isActive( const unsigned& ind ) const ;
///
  unsigned getNumberActive() const ;
///
  unsigned getActiveIndex( const unsigned& ) const ;
/// Transfer derivatives to buffer
  void chainRule( const unsigned&, const unsigned&, const unsigned&, const unsigned&, const double&, const unsigned&, std::vector<double>& buffer );
///
  void copyValues( MultiValue& ) const ;
///
  void copyDerivatives( MultiValue& );
///
  void quotientRule( const unsigned& nder, const unsigned& oder );
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
  hasderiv[nderivatives*ival+jder]=true; hasDerivatives.activate(jder); derivatives[nderivatives*ival+jder] += der;
}

inline
void MultiValue::addTemporyValue( const double& val ) {
  tmpval += val;
}

inline
void MultiValue::addTemporyDerivative( const unsigned& jder, const double& der ) {
  plumed_dbg_assert( jder<nderivatives ); atLeastOneSet=true;
  hasDerivatives.activate(jder); tmpder[jder] += der;
}


inline
void MultiValue::setDerivative( const std::size_t& ival, const std::size_t& jder, const double& der) {
  plumed_dbg_assert( ival<=values.size() && jder<nderivatives ); atLeastOneSet=true;
  hasderiv[nderivatives*ival+jder]=true; hasDerivatives.activate(jder); derivatives[nderivatives*ival+jder]=der;
}


inline
double MultiValue::getDerivative( const std::size_t& ival, const std::size_t& jder ) const {
  plumed_dbg_assert( ival<values.size() && jder<nderivatives );
  // if( !hasderiv[nderivatives*ival+jder] ) return 0.0;
  return derivatives[nderivatives*ival+jder];
}

inline
double MultiValue::getTemporyDerivative( const unsigned& jder ) const {
  plumed_dbg_assert( jder<nderivatives && hasDerivatives.isActive(jder) );
  return tmpder[jder];
}

inline
bool MultiValue::updateComplete() {
  return hasDerivatives.updateComplete();
}

inline
void MultiValue::emptyActiveMembers() {
  hasDerivatives.emptyActiveMembers();
}

inline
void MultiValue::putIndexInActiveArray( const unsigned& ind ) {
  hasDerivatives.putIndexInActiveArray( ind );
}

inline
void MultiValue::updateIndex( const unsigned& ind ) {
  if( hasDerivatives.isActive(ind) ) hasDerivatives.putIndexInActiveArray( ind );
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
bool MultiValue::inVectorCall() const {
  return (nmatrix_cols>0 && vector_call);
}

inline
void MultiValue::clearActiveMembers( const std::size_t& ival ) {
  nactive[ival]=0;
}

inline
void MultiValue::sortActiveList() {
  hasDerivatives.sortActiveList();
}

inline
void MultiValue::completeUpdate() {
  hasDerivatives.completeUpdate();
}

inline
unsigned MultiValue::getNumberActive() const {
  return hasDerivatives.getNumberActive();
}

inline
unsigned MultiValue::getActiveIndex( const unsigned& ind ) const {
  plumed_dbg_assert( ind<hasDerivatives.getNumberActive() );
  return hasDerivatives[ind];
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
void MultiValue::updateDynamicList() {
  if( atLeastOneSet ) hasDerivatives.updateActiveMembers();
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
bool MultiValue::isActive( const unsigned& ind ) const {
  return hasDerivatives.isActive( ind );
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
std::vector<Tensor>& MultiValue::getFirstAtomVirialVector() {
  return tmp_atom_virial;
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
