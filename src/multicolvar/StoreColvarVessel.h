/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#ifndef __PLUMED_multicolvar_StoreColvarVessel_h
#define __PLUMED_multicolvar_StoreColvarVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "vesselbase/StoreDataVessel.h"

namespace PLMD {
namespace multicolvar {

class MultiColvarFunction;

class StoreColvarVessel : public vesselbase::StoreDataVessel {
public:
  static void registerKeywords( Keywords& keys );
  StoreColvarVessel( const vesselbase::VesselOptions& );
  double getValue( const unsigned& );
  virtual std::string description(){ return ""; }
  void chainRuleForComponent( const unsigned& , const unsigned& , const unsigned& , const double& , MultiColvarFunction* );
};

inline
double StoreColvarVessel::getValue( const unsigned& ival ){
  return getComponent( ival, 0 );
}

}
}
#endif
