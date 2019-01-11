/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#ifndef __PLUMED_vesselbase_BridgeVessel_h
#define __PLUMED_vesselbase_BridgeVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "Vessel.h"
#include "core/Value.h"
#include "tools/MultiValue.h"

namespace PLMD {
namespace vesselbase {

/**
\ingroup TOOLBOX
This class allows you to calculate the vessel in one ActionWithVessel.  The user thinks
it is created in a different Action however.  At the moment this is used for region
*/

class BridgeVessel : public Vessel {
private:
  unsigned inum;
  // bool in_normal_calculate;
  std::vector<double> mynumerical_values;
  ActionWithVessel* myOutputAction;
  ActionWithValue* myOutputValues;
  // We create a tempory multivalue here so as to avoid vector resizing
  MultiValue my_tmp_val;
public:
  explicit BridgeVessel( const VesselOptions& );
/// Does this have derivatives
  bool hasDerivatives();
/// Resize the quantities in the vessel
  void resize();
/// Get the action that reads the command in
  ActionWithVessel* getOutputAction();
/// Setup the action we are outputting to
  void setOutputAction( ActionWithVessel* myOutputAction );
/// Apply some force
  bool applyForce( std::vector<double>& forces );
/// Should not be called
  std::string description();
/// Jobs to do before the task list starts
  void prepare();
/// Set the start of the buffer
  void setBufferStart( unsigned& start );
/// This transforms the derivatives using the output value
  MultiValue& transformDerivatives( const unsigned& current, MultiValue& invals, MultiValue& outvals );
/// Actually do the calculation
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_index ) const ;
/// Finish the calculation
  void finish( const std::vector<double>& buffer );
/// Calculate numerical derivatives
  void completeNumericalDerivatives();
/// Set the task flags in the bridged class the same as in the original class
  void copyTaskFlags();
/// Return a tempory multi value - we do this so as to avoid vector resizing
  MultiValue& getTemporyMultiValue();
};

inline
ActionWithVessel* BridgeVessel::getOutputAction() {
  return myOutputAction;
}

}
}
#endif


