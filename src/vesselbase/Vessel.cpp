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
#include "ActionWithVessel.h"
#include "Vessel.h"
#include "tools/Exception.h"
#include "core/Value.h"
#include "tools/Log.h"

namespace PLMD {
namespace vesselbase{

VesselOptions::VesselOptions(const std::string& thisname, const unsigned& nlab, const std::string& params, ActionWithVessel* aa ):
myname(thisname),
numlab(nlab),
action(aa),
parameters(params)
{
} 

Vessel::Vessel( const VesselOptions& da ):
myname(da.myname),
label("unset"),
numlab(da.numlab),
action(da.action),
log((da.action)->log),
comm((da.action)->comm),
serial((da.action)->serial)
{
}

void Vessel::setLabel( const std::string& mylab ){
  plumed_massert( label=="unset", "label has already been set");
  label=mylab;
}

bool Vessel::getLabel( std::string& mylab ) const {
  if( label=="unset" ) return false;
  mylab=label;
  return true;
}

void Vessel::error( const std::string& msg ){
  action->log.printf("ERROR for keyword %s in action %s with label %s : %s \n \n",myname.c_str(), (action->getName()).c_str(), (action->getLabel()).c_str(), msg.c_str() );
  printKeywords();
  plumed_merror("ERROR for keyword " + myname + " in action "  + action->getName() + " with label " + action->getLabel() + " : " + msg );
}

}
}
