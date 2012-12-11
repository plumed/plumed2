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
#ifndef __PLUMED_vessel_base_VesselAccumulator_h
#define __PLUMED_vessel_base_VesselAccumulator_h

#include <string>
#include <cstring>
#include <vector>
#include "VesselValueAccess.h"

namespace PLMD {
namespace vesselbase {

class VesselAccumulator : public VesselValueAccess {
private:
/// The number of buffered values
  unsigned nbuffers;
/// These are pointers to the values in ActionWithValue
  std::vector<Value*> final_values;
protected:
/// Create a value that can be passed between actions
  void addOutput(const std::string& label);
/// Add a value to the buffer
  void addBufferedValue();
/// Get the number of values we are calculating
  unsigned getNumberOfValues() const ;
/// Get pointer to final value
  Value* getPntrToOutput( const unsigned& i );
public:
  VesselAccumulator( const VesselOptions& da );
/// This does the resizing of the buffer
  void resize();
/// This applies all the forces
  bool applyForce( std::vector<double>& forces );
};

inline
Value* VesselAccumulator::getPntrToOutput( const unsigned& iout ){
  plumed_assert( iout<final_values.size() );
  return final_values[iout];
}

inline
unsigned VesselAccumulator::getNumberOfValues() const {
  return final_values.size();
}


}
}
#endif


