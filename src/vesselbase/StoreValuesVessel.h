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
#ifndef __PLUMED_vesselbase_StoreValuesVessel_h
#define __PLUMED_vesselbase_StoreValuesVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "Vessel.h"

namespace PLMD {
namespace vesselbase{

class StoreValuesVessel : public Vessel {
private:
  unsigned bufsize;
  std::vector<unsigned> start;
protected:
/// Are the weights differentiable
  bool diffweight;
/// Get the ith value in the vessel
  double getValue( const unsigned& ) const ;  
/// Get the weight of the ith quantity in the vessel
  double getWeight( const unsigned& ) const ;
/// Add the derivatives from the value
  void addDerivatives( const unsigned& ival, double& pref, Value* value_out );
/// Add the derivatives from the weight of the value
  void addWeightDerivatives( const unsigned& ival, double& pref, Value* value_out );
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  StoreValuesVessel( const VesselOptions& );
/// Return the number of terms
  unsigned getNumberOfTerms(){ return 2; }
/// This does the resizing of the buffer
  void resize();
/// This makes sure all values are stored
  bool calculate();
/// This makes sure things further down the chain are resized
  virtual void local_resizing()=0;
};

inline
double StoreValuesVessel::getValue( const unsigned& ival ) const {
  plumed_dbg_assert( ival<start.size()-1 );
  return getBufferElement( start[ival] );
}

inline
double StoreValuesVessel::getWeight( const unsigned& ival ) const {
  plumed_dbg_assert( ival<start.size()-1 );
  return getBufferElement( bufsize + start[ival] );
}

}
}
#endif
