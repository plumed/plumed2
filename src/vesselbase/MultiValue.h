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
#ifndef __PLUMED_vesselbase_MultiValue_h
#define __PLUMED_vesselbase_MultiValue_h

#include <vector>
#include "tools/Exception.h"
#include "tools/Matrix.h"
#include "tools/DynamicList.h"

namespace PLMD{
namespace vesselbase{

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
/// Set value numbered
  void setValue( const unsigned&,  const double& );
/// Add derivative
  void addDerivative( const unsigned& , const unsigned& , const double& );
/// Return the ith value
  double get( const unsigned& ) const ;
/// Clear all values
  void clearAll();
/// Clear a value
  void clear( const unsigned& );
/// Transfer derivatives to buffer
  void chainRule( const unsigned& , const unsigned& , const unsigned&, const unsigned& , const double& , const unsigned& , std::vector<double>& buffer );
};

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
void MultiValue::addDerivative( const unsigned& ival, const unsigned& jder, const double& der){
  plumed_dbg_assert( ival<=values.size() );
  hasDerivatives.activate(jder); derivatives(ival,jder) += der;
}

}
}
#endif
