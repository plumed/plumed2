/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#ifndef __PLUMED_manyrestraints_ManyRestraintsBase_h
#define __PLUMED_manyrestraints_ManyRestraintsBase_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "vesselbase/ActionWithVessel.h"
#include "vesselbase/ActionWithInputVessel.h"

namespace PLMD {
namespace manyrestraints {

class ManyRestraintsBase :
 public ActionWithValue,
 public ActionPilot,
 public vesselbase::ActionWithVessel,
 public vesselbase::ActionWithInputVessel
{
private:
/// Pointer to underlying action with vessel
  vesselbase::ActionWithVessel* aves;
protected:
/// Get the value of the current cv
  double getValue();
/// Get the weight of the current cv
  double getWeight();
/// Apply the chain rule to calculate the derivatives
  void applyChainRuleForDerivatives( const double& df );
public:
  static void registerKeywords( Keywords& keys );
  ManyRestraintsBase(const ActionOptions&);
  bool isPeriodic(){ return false; }
  unsigned getNumberOfDerivatives();
/// Routines that have to be defined so as not to have problems with virtual methods
  void deactivate_task(){};
/// Don't actually clear the derivatives when this is called from plumed main.  
/// They are calculated inside another action and clearing them would be bad  
  void clearDerivatives(){}
/// Do jobs required before tasks are undertaken
  void doJobsRequiredBeforeTaskList();
// Calculate does nothing
  void calculate(){};
/// Deactivate task now does nothing
  void apply();
  void applyBridgeForces( const std::vector<double>& bb ){ plumed_assert( bb.size()==0 ); }
};

inline
unsigned ManyRestraintsBase::getNumberOfDerivatives(){
  return aves->getNumberOfDerivatives();
}

inline
double ManyRestraintsBase::getValue(){
  return aves->getElementValue(0);
}

inline
double ManyRestraintsBase::getWeight(){
  return aves->getElementValue(1);
}

}
}

#endif
