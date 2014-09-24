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
#include "PathBase.h"
#include "tools/SwitchingFunction.h"

namespace PLMD {
namespace mapping{

void PathBase::registerKeywords( Keywords& keys ){
  Mapping::registerKeywords( keys ); 
  keys.add("compulsory","LAMBDA","the value of the lambda parameter for paths");
  keys.addFlag("NOZPATH",false,"do not calculate the zpath position");
}

PathBase::PathBase(const ActionOptions& ao):
Action(ao),
Mapping(ao)
{
  bool noz; parseFlag("NOZPATH",noz);
  parse("LAMBDA",lambda);

  // Create the list of tasks
  for(unsigned i=0;i<getNumberOfReferencePoints();++i) addTaskToList( i );

  std::string empty="LABEL=zpath";
  if(!noz) addVessel("ZPATH",empty,0);
}

double PathBase::getLambda(){
  return lambda;
}

void PathBase::calculate(){
  // Loop over all frames is now performed by ActionWithVessel
  runAllTasks();
}

void PathBase::performTask(){
  // Calculate the distance from the frame
  double val=calculateDistanceFunction( getCurrentTask(), true );
  // Put the element value in element zero
  setElementValue( 0, val ); setElementValue( 1, 1.0 );
  return;
}

double PathBase::transformHD( const double& dist, double& df ){
  double val = exp( -dist*lambda );
  df = -lambda*val; 
  return val;
}

}
}
