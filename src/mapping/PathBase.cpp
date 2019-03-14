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
#include "PathBase.h"
#include "tools/SwitchingFunction.h"

namespace PLMD {
namespace mapping {

void PathBase::registerKeywords( Keywords& keys ) {
  Mapping::registerKeywords( keys );
  keys.add("compulsory","LAMBDA","0","the value of the lambda parameter for paths");
  keys.addFlag("NOZPATH",false,"do not calculate the zpath position");
}

PathBase::PathBase(const ActionOptions& ao):
  Action(ao),
  Mapping(ao)
{
  setLowMemOption( true );
  weightHasDerivatives=true;
  bool noz; parseFlag("NOZPATH",noz);
  parse("LAMBDA",lambda);

  // Create the list of tasks
  for(unsigned i=0; i<getNumberOfReferencePoints(); ++i) addTaskToList( i );
  // And activate them all
  deactivateAllTasks();
  for(unsigned i=0; i<getFullNumberOfTasks(); ++i) taskFlags[i]=1;
  lockContributors();

  std::string empty="LABEL=zpath";
  if(!noz) {
    if( lambda==0 ) error("you must set LAMDBA value in order to calculate ZPATH coordinate.  Use LAMBDA/NOZPATH keyword");
    addVessel("ZPATH",empty,0);
  }
}

double PathBase::getLambda() {
  return lambda;
}

void PathBase::calculate() {
  // Loop over all frames is now performed by ActionWithVessel
  runAllTasks();
}

void PathBase::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  // This builds a pack to hold the derivatives
  ReferenceValuePack mypack( getNumberOfArguments(), getNumberOfAtoms(), myvals );
  finishPackSetup( current, mypack );
  // Calculate the distance from the frame
  double val=calculateDistanceFunction( current, mypack, true );
  // Put the element value in element zero
  myvals.setValue( 0, val ); myvals.setValue( 1, 1.0 );
  return;
}

double PathBase::transformHD( const double& dist, double& df ) const {
  if( lambda==0 ) { df=1; return dist; }
  double val = exp( -dist*lambda );
  df = -lambda*val;
  return val;
}

}
}
