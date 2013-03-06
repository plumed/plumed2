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
#ifndef __PLUMED_vesselbase_BridgeVessel_h
#define __PLUMED_vesselbase_BridgeVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "Vessel.h"
#include "core/Value.h"

namespace PLMD {
namespace vesselbase {

/**
\ingroup TOOLBOX
This class allows you to calculate the vessel in one ActionWithVessel.  The user thinks
it is created in a different Action however.  At the moment this is used for region 
*/

class BridgeVessel : public Vessel {
private:
  std::vector<double> forces;
  unsigned inum;
  std::vector<double> mynumerical_values;
  ActionWithVessel* myOutputAction;
  ActionWithValue* myOutputValues;
public:
  BridgeVessel( const VesselOptions& );
/// Resize the quantities in the vessel
  void resize();
/// Setup the action we are outputting to
  void setOutputAction( ActionWithVessel* myOutputAction );
/// Apply some force 
  bool applyForce( std::vector<double>& forces );
/// Should not be called
  std::string description();
/// Want to get rid of this function
  unsigned getNumberOfTerms(){ return 2; }
/// Jobs to do before the task list starts
  void prepare();
/// Actually do the calculation
  bool calculate();
/// Finish the calculation
  void finish();
/// Calculate numerical derivatives
  void completeNumericalDerivatives();
};

}
}
#endif


