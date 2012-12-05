/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_VesselValueAccess_h
#define __PLUMED_VesselValueAccess_h

#include <string>
#include <cstring>
#include <vector>
#include "Vessel.h"
#include "core/Value.h"

namespace PLMD {

class VesselValueAccess : public Vessel {
private:
/// The value we take from the action
  Value myvalue;
protected:
/// The start for each of the values in the system
  std::vector<unsigned> value_starts;
/// Set the number of values
  void setNumberOfValues( const unsigned& );
/// Set the sizes of all the values
  void setValueSizes( const std::vector<unsigned>& );
/// Get the ith value from the buffer
  void getValue( const unsigned& , Value& ) const ;
/// Get the value of the ith value in the buffer
  double getValue( const unsigned& ) const ;
public:
/// Constructor
  VesselValueAccess( const VesselOptions& );
};

inline
void VesselValueAccess::getValue( const unsigned& icv, Value& val ) const {
   plumed_assert( icv<value_starts.size()-1 );
   unsigned nder=(value_starts[icv+1]-value_starts[icv]-1);
   if( val.getNumberOfDerivatives()!=nder ) val.resizeDerivatives( nder );
   val.clearDerivatives();
   unsigned ider=value_starts[icv]; val.set( getBufferElement(ider) ); ider++;
   for(unsigned i=0;i<nder;++i){ val.addDerivative( i, getBufferElement(ider) ); ider++; }
}

inline
double VesselValueAccess::getValue( const unsigned& icv ) const {
   plumed_assert( icv<value_starts.size()-1 );
   return getBufferElement( value_starts[icv] );
}

}  
#endif
