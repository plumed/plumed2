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
public:
  static void registerKeywords( Keywords& keys );
  explicit ManyRestraintsBase(const ActionOptions&);
  bool isPeriodic() override { return false; }
  unsigned getNumberOfDerivatives() override;
/// Routines that have to be defined so as not to have problems with virtual methods
  void deactivate_task( const unsigned & task_index ) {};
/// Don't actually clear the derivatives when this is called from plumed main.
/// They are calculated inside another action and clearing them would be bad
  void clearDerivatives() override {}
/// Do jobs required before tasks are undertaken
  void doJobsRequiredBeforeTaskList() override;
/// This actually does the calculation
  void transformBridgedDerivatives( const unsigned& current, MultiValue& invals, MultiValue& outvals ) const override;
/// Calculate the potential
  virtual double calcPotential( const double& val, double& df ) const=0;
// Calculate does nothing
  void calculate() override {};
/// This should never be called
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const override { plumed_error(); }
/// Deactivate task now does nothing
  void apply() override;
  void applyBridgeForces( const std::vector<double>& bb ) override { plumed_assert( bb.size()==0 ); }
};

inline
unsigned ManyRestraintsBase::getNumberOfDerivatives() {
  return aves->getNumberOfDerivatives();
}

}
}

#endif
