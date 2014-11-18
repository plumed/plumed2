/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#ifndef __PLUMED_tools_MultiValue_h
#define __PLUMED_tools_MultiValue_h

#include <vector>
#include "Exception.h"
#include "Matrix.h"
#include "DynamicList.h"

namespace PLMD{

class MultiValue {
private:
/// Used to ensure rapid accumulation of derivatives
  DynamicList<unsigned> hasDerivatives;
/// Values of quantities
  std::vector<double> values;
/// Derivatives
  Matrix<double> derivatives;
public:
  MultiValue( const unsigned& , const unsigned& );
  void resize( const unsigned& , const unsigned& );
/// Get the number of values in the stash
  unsigned getNumberOfValues() const ; 
/// Get the number of derivatives in the stash
  unsigned getNumberOfDerivatives() const ;
/// Set value numbered
  void setValue( const unsigned&,  const double& );
/// Add value numbered
  void addValue( const unsigned&,  const double& );
/// Add derivative
  void addDerivative( const unsigned& , const unsigned& , const double& );
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
/// Functions for accessing active list
  bool updateComplete();
  void emptyActiveMembers();
  void updateIndex( const unsigned & );
  void sortActiveList();
  void updateDynamicList();
///
  unsigned getNumberActive() const ;
///
  unsigned getActiveIndex( const unsigned& ) const ;
/// Transfer derivatives to buffer
  void chainRule( const unsigned& , const unsigned& , const unsigned&, const unsigned& , const double& , const unsigned& , std::vector<double>& buffer );
///
  void copyValues( MultiValue& ) const ;
///
  void copyDerivatives( MultiValue& );
///
  void quotientRule( const unsigned& nder, const unsigned& dder, const unsigned& oder );
};

inline
unsigned MultiValue::getNumberOfValues() const {
  return values.size();
}

inline
unsigned MultiValue::getNumberOfDerivatives() const {
  return derivatives.ncols();
}

inline
double MultiValue::get( const unsigned& ival ) const {
  plumed_dbg_assert( ival<=values.size() );
  return values[ival];
}

inline
void MultiValue::setValue( const unsigned& ival,  const double& val){
  plumed_dbg_assert( ival<=values.size() );
  values[ival]=val;
}

inline
void MultiValue::addValue( const unsigned& ival,  const double& val){
  plumed_dbg_assert( ival<=values.size() );
  values[ival]+=val;
}

inline
void MultiValue::addDerivative( const unsigned& ival, const unsigned& jder, const double& der){
  plumed_dbg_assert( ival<=values.size() );
  hasDerivatives.activate(jder); derivatives(ival,jder) += der;
}

inline
void MultiValue::setDerivative( const unsigned& ival, const unsigned& jder, const double& der){
  plumed_dbg_assert( ival<=values.size() );
  hasDerivatives.activate(jder); derivatives(ival,jder)=der;
}


inline
double MultiValue::getDerivative( const unsigned& ival, const unsigned& jder ) const {
  plumed_dbg_assert( jder<derivatives.ncols() && hasDerivatives.isActive(jder) );
  return derivatives(ival,jder);
}

inline
bool MultiValue::updateComplete(){
  return hasDerivatives.updateComplete();
}

inline
void MultiValue::emptyActiveMembers(){
  hasDerivatives.emptyActiveMembers();
}

inline
void MultiValue::updateIndex( const unsigned& ind ){
  hasDerivatives.updateIndex( ind );
}

inline
void MultiValue::sortActiveList(){
  hasDerivatives.sortActiveList();
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
void MultiValue::updateDynamicList(){
  hasDerivatives.updateActiveMembers();
}

}
#endif
