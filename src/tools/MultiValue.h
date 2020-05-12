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

#include <vector>
#include "Exception.h"
#include "DynamicList.h"

namespace PLMD {

class MultiValue {
private:
/// Used to ensure rapid accumulation of derivatives
  DynamicList<unsigned> hasDerivatives;
/// Values of quantities
  std::vector<double> values;
/// Number of derivatives per value
  unsigned nderivatives;
/// Derivatives
  std::vector<double> derivatives;
/// Tempory value
  double tmpval;
/// Tempory vector of derivatives (used for calculating quotients
  std::vector<double> tmpder;
/// Logical to check if any derivatives were set
  bool atLeastOneSet;
/// This is a fudge to save on vector resizing in MultiColvar
  std::vector<unsigned> indices, sort_indices;
  std::vector<Vector> tmp_atoms;
public:
  MultiValue( const unsigned&, const unsigned& );
  void resize( const unsigned&, const unsigned& );
///
  std::vector<unsigned>& getIndices();
  std::vector<unsigned>& getSortIndices();
  std::vector<Vector>& getAtomVector();
/// Get the number of values in the stash
  unsigned getNumberOfValues() const ;
/// Get the number of derivatives in the stash
  unsigned getNumberOfDerivatives() const ;
/// Set value numbered
  void setValue( const unsigned&,  const double& );
/// Add value numbered
  void addValue( const unsigned&,  const double& );
/// Add derivative
  void addDerivative( const unsigned&, const unsigned&, const double& );
/// Add to the tempory value
  void addTemporyValue( const double& val );
/// Add tempory derivatives - this is used for calculating quotients
  void addTemporyDerivative( const unsigned& jder, const double& der );
/// Set the value of the derivative
  void setDerivative( const unsigned& ival, const unsigned& jder, const double& der);
/// Return the ith value
  double get( const unsigned& ) const ;
/// Return a derivative value
  double getDerivative( const unsigned&, const unsigned& ) const ;
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
double MultiValue::get( const unsigned& ival ) const {
  plumed_dbg_assert( ival<=values.size() );
  return values[ival];
}

inline
void MultiValue::setValue( const unsigned& ival,  const double& val) {
  plumed_dbg_assert( ival<=values.size() );
  values[ival]=val;
}

inline
void MultiValue::addValue( const unsigned& ival,  const double& val) {
  plumed_dbg_assert( ival<=values.size() );
  values[ival]+=val;
}

inline
void MultiValue::addDerivative( const unsigned& ival, const unsigned& jder, const double& der) {
  plumed_dbg_assert( ival<=values.size() && jder<nderivatives ); atLeastOneSet=true;
  hasDerivatives.activate(jder); derivatives[nderivatives*ival+jder] += der;
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
void MultiValue::setDerivative( const unsigned& ival, const unsigned& jder, const double& der) {
  plumed_dbg_assert( ival<=values.size() && jder<nderivatives ); atLeastOneSet=true;
  hasDerivatives.activate(jder); derivatives[nderivatives*ival+jder]=der;
}


inline
double MultiValue::getDerivative( const unsigned& ival, const unsigned& jder ) const {
  plumed_dbg_assert( jder<nderivatives && hasDerivatives.isActive(jder) );
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

}
#endif
