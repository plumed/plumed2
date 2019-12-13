/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#ifndef __PLUMED_vesselbase_ActionWithInputVessel_h
#define __PLUMED_vesselbase_ActionWithInputVessel_h

#include "core/Action.h"
#include "Vessel.h"
#include <vector>

namespace PLMD {
namespace vesselbase {

/**
\ingroup MULTIINHERIT
*/

class ActionWithInputVessel : public virtual Action {
private:
  Vessel* arguments;
  BridgeVessel* myBridgeVessel;
protected:
/// What type of arguments are we reading in
  void readArgument( const std::string& type );
/// Return a pointer to specific argument
  Vessel* getPntrToArgument();
/// Add forces to arguments (used in apply)
  void addForcesOnArguments( const std::vector<double>& forces );
public:
/// Registers the list of keywords
  static void registerKeywords( Keywords& keys );
  explicit ActionWithInputVessel(const ActionOptions&);
  virtual ~ActionWithInputVessel() {}
/// Calculate the numerical derivatives
/// N.B. only pass an ActionWithValue to this routine if you know exactly what you
/// are doing.  The default will be correct for the vast majority of cases
  virtual void calculateNumericalDerivatives( ActionWithValue* a=NULL );
/// Apply forces from the bridge
  void applyBridgeForces( const std::vector<double>& bb );
/// Apply forces from the bridge
  virtual void addBridgeForces( const std::vector<double>& bb ) {}
};

inline
Vessel* ActionWithInputVessel::getPntrToArgument() {
  return arguments;
}

}
}

#endif
