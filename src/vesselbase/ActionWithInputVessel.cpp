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
#include "ActionWithInputVessel.h"
#include "StoreDataVessel.h"
#include "BridgeVessel.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace vesselbase {

void ActionWithInputVessel::registerKeywords(Keywords& keys) {
  keys.add("compulsory","DATA","certain actions in plumed work by calculating a list of variables and summing over them. "
           "This particular action can be used to calculate functions of these base variables or prints "
           "them to a file. This keyword thus takes the label of one of those such variables as input.");
  keys.reserve("compulsory","GRID","the action that creates the input grid you would like to use");
}

ActionWithInputVessel::ActionWithInputVessel(const ActionOptions&ao):
  Action(ao),
  arguments(NULL),
  myBridgeVessel(NULL)
{
}

void ActionWithInputVessel::readArgument( const std::string& type ) {
  std::string mlab;
  if( keywords.exists("DATA") && type!="grid" ) parse("DATA",mlab);
  ActionWithVessel* mves= plumed.getActionSet().selectWithLabel<ActionWithVessel*>(mlab);
  if(!mves) error("action labelled " +  mlab + " does not exist or does not have vessels");
  addDependency(mves);

  ActionWithValue* aval=dynamic_cast<ActionWithValue*>( this );
  if(aval) {
    if( aval->checkNumericalDerivatives() ) {
      ActionWithValue* aval2=dynamic_cast<ActionWithValue*>( mves );
      plumed_assert( aval2 ); aval2->useNumericalDerivatives();
    }
  }

  if( type=="bridge" ) {
    ActionWithVessel* aves=dynamic_cast<ActionWithVessel*>( this );
    plumed_assert(aves); myBridgeVessel = mves->addBridgingVessel( aves );
    arguments = dynamic_cast<Vessel*>( myBridgeVessel );
  } else  if( type=="store" ) {
    arguments = dynamic_cast<Vessel*>( mves->buildDataStashes( NULL ) );
  } else {
    plumed_error();
  }
}

void ActionWithInputVessel::calculateNumericalDerivatives( ActionWithValue* a ) {
  if(!a) {
    a=dynamic_cast<ActionWithValue*>(this);
    plumed_massert(a,"cannot compute numerical derivatives for an action without values");
  }
  if( myBridgeVessel ) {
    myBridgeVessel->completeNumericalDerivatives();
  } else {
    error("numerical derivatives are not implemented");
  }
}

void ActionWithInputVessel::applyBridgeForces( const std::vector<double>& bb ) {
  plumed_dbg_assert( myBridgeVessel ); addBridgeForces( bb );
}

}
}
