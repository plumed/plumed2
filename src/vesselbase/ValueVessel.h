/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#ifndef __PLUMED_vesselbase_ValueVessel_h
#define __PLUMED_vesselbase_ValueVessel_h

#include <string>
#include <cstring>
#include <vector>

#include "Vessel.h"
#include "core/Value.h"

namespace PLMD {
namespace vesselbase{

class ValueVessel : public Vessel {
private: 
  bool no_output_value;
  Value* final_value;
protected:
/// The component that is being averaged/accumulated whatever
  unsigned mycomp;
/// Set the final value
  void setOutputValue( const double& val );
public:
  static void registerKeywords( Keywords& keys );
  explicit ValueVessel( const VesselOptions& da );
  ~ValueVessel(); 
  std::string description();
  virtual std::string value_descriptor()=0;
  bool applyForce( std::vector<double>& forces );
  double getOutputValue() const ;
  Value* getFinalValue() const ;
};

inline
Value* ValueVessel::getFinalValue() const {
  return final_value;
}

inline
double ValueVessel::getOutputValue() const {
  return final_value->get();
}

inline
void ValueVessel::setOutputValue( const double& val ){
  final_value->set( val );
}

}
}
#endif
