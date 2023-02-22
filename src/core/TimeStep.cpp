/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "ActionToPutData.h"
#include "ActionRegister.h"

namespace PLMD {

class TimeStep : public ActionToPutData {
private: 
  double timestep;
public:
  explicit TimeStep(const ActionOptions&);
// active methods:
  static void registerKeywords( Keywords& keys );
/// Set the value pointer
  bool setValuePointer( const std::string& name, void* val ) override ;
  void transferFixedValue( const double& unit ) override ;
};


PLUMED_REGISTER_ACTION(TimeStep,"TIMESTEP")

TimeStep::TimeStep(const ActionOptions&ao):
  Action(ao),
  ActionToPutData(ao),
  timestep(0)
{
  std::vector<unsigned> shape; addValue( shape ); 
  setNotPeriodic(); setUnit( "time", "default" );
}

void TimeStep::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); 
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
}

bool TimeStep::setValuePointer( const std::string& name, void* val ) {
  if( name!=getLabel() ) return false;
  wasset=true; plumed_massert( dataCanBeSet, "set " + getLabel() + " cannot be set at this time"); 
  timestep = MD2double(val);  
// The following is to avoid extra digits in case the MD code uses floats
// e.g.: float f=0.002 when converted to double becomes 0.002000000094995
// To avoid this, we keep only up to 6 significant digits after first one
  double magnitude=std::pow(10,std::floor(std::log10(timestep)));
  timestep=std::floor(timestep/magnitude*1e6)/1e6*magnitude;
  return true; 
}

void TimeStep::transferFixedValue( const double& unit ) {
   plumed_assert( fixed ); if( !wasset ) return;
// We have to transfer the value here as we store the data as a double in PLUMED in order to 
// overcome all the stuff about extra digits in the above comment
   getPntrToValue()->set( unit*timestep );
}  


}



