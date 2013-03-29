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
#include "MultiColvarFunction.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "MultiColvar.h"

namespace PLMD {
namespace multicolvar { 

void MultiColvarFunction::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.add("compulsory","ARG","the label of the action that calculates the vectors we are interested in");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("optional","NL_STRIDE","the frequency with which the neighbor list should be updated. Between neighbour list update steps all quantities "
                                  "that contributed less than TOL at the previous neighbor list update step are ignored.");
  vesselbase::ActionWithVessel::registerKeywords( keys );
}

MultiColvarFunction::MultiColvarFunction(const ActionOptions& ao):
Action(ao),
ActionWithValue(ao),
ActionWithVessel(ao),
usepbc(false),
updateFreq(0),
lastUpdate(0)
{
  weightHasDerivatives=true; // In most cases the weight will have derivatives

  std::string mlab; parse("ARG",mlab);
  mycolv = plumed.getActionSet().selectWithLabel<multicolvar::MultiColvar*>(mlab); 
  if(!mycolv) error("action labeled " + mlab + " does not exist or is not a multicolvar");
  addDependency(mycolv);

  // Retrieve the central atoms
  catoms = mycolv->getCentralAtoms();

  if( keywords.exists("NOPBC") ){
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  }
  if( keywords.exists("NL_STRIDE") ) parse("NL_STRIDE",updateFreq);
  if(updateFreq>0) log.printf("  Updating contributors every %d steps.\n",updateFreq);
  else log.printf("  Updating contributors every step.\n");
}

void MultiColvarFunction::addColvar( const std::vector<unsigned>& newatoms ){
  if( colvar_atoms.size()>0 ) plumed_assert( colvar_atoms[0].fullSize()==newatoms.size() );
  DynamicList<unsigned> newlist;
  for(unsigned i=0;i<newatoms.size();++i) newlist.addIndexToList( newatoms[i] );
  taskList.addIndexToList( colvar_atoms.size() );
  colvar_atoms.push_back( newlist );
}

void MultiColvarFunction::completeSetup(){
  taskList.activateAll();
  for(unsigned i=0;i<colvar_atoms.size();++i) colvar_atoms[i].activateAll();
}

void MultiColvarFunction::prepare(){
  bool updatetime=false;
  if( contributorsAreUnlocked ){
      taskList.mpi_gatherActiveMembers( comm );
      mpi_gatherActiveMembers( comm, colvar_atoms ); 
      lockContributors(); updatetime=true;
  }
  if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
      taskList.activateAll(); 
      for(unsigned i=0;i<taskList.getNumberActive();++i) colvar_atoms[i].activateAll();
      unlockContributors(); updatetime=true; lastUpdate=getStep();
  }
  if(updatetime) resizeFunctions();
}

Vector MultiColvarFunction::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc){ return mycolv->pbcDistance( vec1, vec2 ); }
  else{ return delta( vec1, vec2 ); }
}

void MultiColvarFunction::calculate(){
  if( checkNumericalDerivatives() ){
     // Clear any derivatives in base colvar that were 
     // accumulated from previous calculations
     mycolv->clearDerivatives(); 
     // And recalculate
     mycolv->calculate();
  }
  runAllTasks();
}

void MultiColvarFunction::calculateNumericalDerivatives( ActionWithValue* a ){
  mycolv->calculateNumericalDerivatives( this );
}

void MultiColvarFunction::performTask( const unsigned& j ){
  // Compute a weight for this quantity
  double weight=calculateWeight();
  setElementValue( 1, weight );
  if( weight<getTolerance() ) return;

  // Now compute the quantity of interest
  double val=compute(); 

  // And set everything ready for vessels
  setElementValue( 0, weight*val );
  if( weightHasDerivatives ){
      unsigned nder=getNumberOfDerivatives();
      for(unsigned i=mycolv->getFirstDerivativeToMerge();i<mycolv->getNumberOfDerivatives();i=mycolv->getNextDerivativeToMerge(i)){
          setElementDerivative( i, weight*getElementDerivative(i) + val*getElementDerivative(nder+i) );
      }
  }
  return;
}

void MultiColvarFunction::apply(){
  std::vector<double> forces( getNumberOfDerivatives(), 0.0 );
  for(int i=0;i<getNumberOfVessels();++i){
     if( (getPntrToVessel(i)->applyForce( forces )) ){
          catoms->addForces( forces );
     } 
  }
}

}
}

