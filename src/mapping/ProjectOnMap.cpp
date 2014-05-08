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
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "tools/ConjugateGradient.h"
#include "vesselbase/ActionWithInputVessel.h"
#include "Mapping.h"

namespace PLMD {
namespace mapping {

class ProjectOnMap :
  public ActionWithValue,
  public vesselbase::ActionWithInputVessel
{
private:
  Mapping* mymap;
  double tolerance;
  std::vector<double> mypoint;
public:
  static void registerKeywords( Keywords& keys );
  ProjectOnMap(const ActionOptions& ao);
  unsigned getNumberOfDerivatives(){ return 0; }
  void calculate();
  void apply(){}
};

PLUMED_REGISTER_ACTION(ProjectOnMap,"FINDPROJECTION")

void ProjectOnMap::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  vesselbase::ActionWithInputVessel::registerKeywords( keys );
  keys.remove("DATA"); keys.use("FUNC");
//  keys.add("compulsory","ARG","The name of the action that describes the mapping from the low to the high dimensional space");
  keys.add("compulsory","TOL","1E-6","The tolerance for the conjugate gradient optimization algorithm");
}

ProjectOnMap::ProjectOnMap(const ActionOptions& ao):
Action(ao),
ActionWithValue(ao),
ActionWithInputVessel(ao),
mymap(NULL)
{
  // Find the mapping object
  readArgument( "func" ); plumed_assert( getDependencies().size()==1 );
  mymap=dynamic_cast<Mapping*>( getDependencies()[0] );
//  std::string mylab; parse("ARG",mylab);
//  mymap=plumed.getActionSet().selectWithLabel<Mapping*>(mylab);
//  if(!mymap) error(mylab + " mapping does not exist");
//  addDependency(mymap);

  // Read in the tolerance for the CG minimisation
  parse("TOL",tolerance);

  // Resise the projection coordinates
  mypoint.resize( mymap->getNumberOfProperties() );
  // Create components
  for(unsigned i=0;i<mypoint.size();++i){
     addComponent( mymap->getPropertyName(i) + "_proj" );
     componentIsNotPeriodic( mymap->getPropertyName(i) + "_proj" );
  }
}

void ProjectOnMap::calculate(){
  // Find the start point (projection of point closest to the current position)
  mymap->findClosestPoint( mypoint );

  // Now do the minimisation
  ConjugateGradient<Mapping> myminimiser( mymap );
  myminimiser.minimise( tolerance, mypoint, &Mapping::calculateStress );

  // And set the values
  for(unsigned i=0;i<mypoint.size();++i) getPntrToComponent(i)->set( mypoint[i] );
}

}
}

