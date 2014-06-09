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
#ifndef __PLUMED_multicolvar_BridgedMultiColvarFunction_h
#define __PLUMED_multicolvar_BridgedMultiColvarFunction_h

#include "core/ActionAtomistic.h"
#include "tools/HistogramBead.h"
#include "tools/Pbc.h"
#include "core/ActionWithValue.h"
#include "vesselbase/ActionWithVessel.h"
#include "vesselbase/ActionWithInputVessel.h"
#include "vesselbase/BridgeVessel.h"
#include "MultiColvarBase.h"

namespace PLMD {
namespace multicolvar {

class BridgedMultiColvarFunction :
  public ActionAtomistic,
  public ActionWithValue,
  public vesselbase::ActionWithVessel,
  public vesselbase::ActionWithInputVessel
  {
friend class MultiColvarBase;   
private:
/// This is used for storing positions properly
  Vector tmp_p;
/// The action that is calculating the colvars of interest
  MultiColvarBase* mycolv;
/// The vessel that bridges
  vesselbase::BridgeVessel* myBridgeVessel;
/// Everything for controlling the updating of neighbor lists
  bool firsttime;
  int updateFreq;
/// Fast merging of derivatives (automatic skips of zero contributions)
  DynamicList<unsigned> activeAtoms;
protected:
/// Get a pointer to the base multicolvar
  MultiColvarBase* getPntrToMultiColvar() const ;
/// Deactivate all the atoms in the list
  void deactivateAllAtoms();
/// Activate the nth atom in the list
  void setAtomActive( const unsigned& n );
/// And complete the list of active atoms
  void updateActiveAtoms();
public:
  static void registerKeywords( Keywords& keys );
  BridgedMultiColvarFunction(const ActionOptions&);
/// Don't actually clear the derivatives when this is called from plumed main.  
/// They are calculated inside another action and clearing them would be bad  
  void clearDerivatives(){}
// This is used during neighbor list update step
  void finishTaskListUpdate();
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Turn on the derivatives
  void turnOnDerivatives();
/// Is the output quantity periodic
  bool isPeriodic();
/// Jobs to be done when the action is activated
  void prepare();
/// Routines that have to be defined so as not to have problems with virtual methods 
  void deactivate_task();
  void calculate(){}
/// We need our own calculate numerical derivatives here
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  void apply(){};
/// These routines replace the virtual routines in ActionWithVessel for 
/// code optimization
  void mergeDerivatives( const unsigned& ider, const double& df );
  void clearDerivativesAfterTask( const unsigned& ider );
};

inline
MultiColvarBase* BridgedMultiColvarFunction::getPntrToMultiColvar() const {
  return mycolv;
}

inline
unsigned BridgedMultiColvarFunction::getNumberOfDerivatives(){
  return mycolv->getNumberOfDerivatives() + 3*getNumberOfAtoms();
}

inline
void BridgedMultiColvarFunction::deactivateAllAtoms(){
  activeAtoms.deactivateAll();
}

inline
void BridgedMultiColvarFunction::setAtomActive( const unsigned& n ){
  activeAtoms.activate(n);
}

inline
void BridgedMultiColvarFunction::updateActiveAtoms(){
  activeAtoms.updateActiveMembers();
}

}
}
#endif
