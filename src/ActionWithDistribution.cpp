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
#include "ActionWithDistribution.h"
#include "MultiColvar.h"

using namespace std;
using namespace PLMD;

void ActionWithDistribution::registerKeywords(Keywords& keys){
  keys.add("optional","TOL","when accumulating sums quantities that contribute less than this will be ignored.");
  keys.add("optional","NL_STRIDE","the frequency with which the neighbor list should be updated. Between neighbour list update steps all quantities that contributed less than TOL at the previous neighbor list update step are ignored.");
}

void ActionWithDistribution::autoParallelize(Keywords& keys){
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize over collective variables");
}

ActionWithDistribution::ActionWithDistribution(const ActionOptions&ao):
  Action(ao),
  read(false),
  all_values(true),
  serial(false),
  updateFreq(0),
  lastUpdate(0),
  reduceAtNextStep(false),
  tolerance(0),
  myfield(NULL)
{
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  else serial=true;
  if(serial)log.printf("  doing calculation in serial\n");
  tolerance=epsilon; 
  if( keywords.exists("TOL") ) parse("TOL",tolerance);
  if( keywords.exists("NL_STRIDE") ) parse("NL_STRIDE",updateFreq);
  if(updateFreq>0){
    log.printf("  Updating contributors every %d steps. Ignoring contributions less than %lf\n",updateFreq,tolerance);
  } else {
    log.printf("  Updating contributors every step.");
    if( tolerance>0) log.printf(" Ignoring contributions less than %lf\n",tolerance);
    else log.printf("\n");
  }
}

ActionWithDistribution::~ActionWithDistribution(){
  delete myfield;
  for(unsigned i=0;i<functions.size();++i) delete functions[i]; 
}

void ActionWithDistribution::addField( std::string key, Field* ff ){
  plumed_assert( key.length()!=0 );
  std::string fieldin; parse(key,fieldin);
  if( fieldin.length()==0 ) return ;
  all_values=false;  std::string freport;
  myfield=ff; myfield->read( fieldin, getNumberOfFunctionsInDistribution(), freport );
  if( !myfield->check() ){
     log.printf("ERROR for keyword FIELD in action %s with label %s : %s \n \n", getName().c_str(), getLabel().c_str(), ( myfield->errorMessage() ).c_str() );
     myfield->printKeywords( log );
     plumed_merror("ERROR for keyword FIELD in action "  + getName() + " with label " + getLabel() + " : " + myfield->errorMessage() );
     exit(1);
  }
  log.printf("  %s\n", freport.c_str() );
  derivedFieldSetup( myfield->get_sigma() );
}

void ActionWithDistribution::addDistributionFunction( std::string name, DistributionFunction* fun ){
  if(all_values) all_values=false;  // Possibly will add functionality to delete all values here

  // Some sanity checks
  gradient* gfun=dynamic_cast<gradient*>(fun);
  cvdens* cvfun=dynamic_cast<cvdens*>(fun);
  if( gfun || cvfun ){
      MultiColvar* mcheck=dynamic_cast<MultiColvar*>(this);
      if(!mcheck) plumed_massert(mcheck,"cannot do gradient or cvdens if this is not a MultiColvar");
  } 

  // Check function is good 
  if( !fun->check() ){
     log.printf("ERROR for keyword %s in action %s with label %s : %s \n \n",name.c_str(), getName().c_str(), getLabel().c_str(), ( fun->errorMessage() ).c_str() );
     fun->printKeywords( log ); 
     plumed_merror("ERROR for keyword " + name + " in action "  + getName() + " with label " + getLabel() + " : " + fun->errorMessage() );
     exit(1);
  }

  // Add a value
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  a->addComponentWithDerivatives( fun->getLabel() );
  a->componentIsNotPeriodic( fun->getLabel() );
  unsigned fno=a->getNumberOfComponents()-1;
  final_values.push_back( a->copyOutput( fno ) );

  // Add the function   
  plumed_massert( fno==functions.size(), "Number of functions does not match number of values" );
  functions.push_back( fun );
  log.printf("  value %s contains %s\n",( (a->copyOutput( fno ))->getName() ).c_str(),( functions[fno]->message() ).c_str() );
}

void ActionWithDistribution::requestDistribution(){
  read=true; 
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  plumed_massert(a,"can only do distribution on ActionsWithValue");

  if(all_values){
      error("No function has been specified");
  } 
  plumed_massert( functions.size()==final_values.size(), "number of functions does not match number of values" );
  // This sets up the dynamic list that holds what we are calculating
  for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){ members.addIndexToList(i); }
  members.activateAll(); members.updateActiveMembers();
  // We prepare the first step as if we are doing a neighbor list update
  prepareForNeighborListUpdate(); reduceAtNextStep=true;
}

void ActionWithDistribution::prepare(){
 if(reduceAtNextStep){
    completeNeighborListUpdate();
    // Setup the functions by declaring enough space to hold the derivatives
    ActionWithValue* av=dynamic_cast<ActionWithValue*>(this);
    for(unsigned i=0;i<functions.size();++i) functions[i]->setNumberOfDerivatives( (av->getPntrToComponent(i))->getNumberOfDerivatives() );
    // Setup the buffers for mpi gather
    if( myfield ){
      std::vector<unsigned> cv_sizes( getNumberOfFunctionsInDistribution() ); 
      for(unsigned i=0;i<cv_sizes.size();++i){ cv_sizes[i]=getThisFunctionsNumberOfDerivatives(i); }
      myfield->resizeBaseQuantityBuffers( cv_sizes ); 
      myfield->resizeDerivatives( getNumberOfFieldDerivatives() );
    } 
    if( functions.size()!=0 ){
      unsigned bufsize=0;
      for(unsigned i=0;i<functions.size();++i) bufsize+=functions[i]->requiredBufferSpace();
      buffer.resize( bufsize ); 
    }
    reduceAtNextStep=false;
 }
 if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
    members.activateAll();
    members.updateActiveMembers();
    prepareForNeighborListUpdate();
    // Setup the functions by declaring enough space to hold the derivatives
    ActionWithValue* av=dynamic_cast<ActionWithValue*>(this);
    for(unsigned i=0;i<functions.size();++i) functions[i]->setNumberOfDerivatives( (av->getPntrToComponent(i))->getNumberOfDerivatives() );
    // Setup the buffers for mpi gather
    if( myfield ){
      std::vector<unsigned> cv_sizes( getNumberOfFunctionsInDistribution() ); unsigned kk;
      for(unsigned i=0;i<cv_sizes.size();++i){ cv_sizes[i]=getThisFunctionsNumberOfDerivatives(i); }
      myfield->resizeBaseQuantityBuffers( cv_sizes ); 
      myfield->resizeDerivatives( getNumberOfFieldDerivatives() );
    } 
    if( functions.size()!=0 ){
      unsigned bufsize=0;
      for(unsigned i=0;i<functions.size();++i) bufsize+=functions[i]->requiredBufferSpace();
      buffer.resize( bufsize ); 
    }
    reduceAtNextStep=true;
    lastUpdate=getStep();
 }
}

void ActionWithDistribution::calculate(){
  plumed_massert( read, "you must have a call to requestDistribution somewhere" );
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){ stride=1; rank=0; }

  std::vector<Value> aux;
  // Reset everything
  for(unsigned j=0;j<functions.size();++j) functions[j]->reset();
  if( myfield ) myfield->clear();
  // Create a value to store stuff in 
  Value* tmpvalue=new Value();

  unsigned kk; bool keep;
  for(unsigned i=rank;i<members.getNumberActive();i+=stride){
      // Retrieve the function we are calculating from the dynamic list
      kk=members[i]; 
      // Make sure we have enough derivatives in this value
      unsigned nder=getThisFunctionsNumberOfDerivatives(kk);
      if( tmpvalue->getNumberOfDerivatives()!=nder ) tmpvalue->resizeDerivatives( nder );
      // Retrieve the periodicity of this value
      if( isPeriodic(kk) ){ 
         double min, max; retrieveDomain( kk, min, max );
         tmpvalue->setDomain( min, max ); 
      } else { tmpvalue->setNotPeriodic(); }
 
      // Calculate the value of this particular function 
      calculateThisFunction( kk, tmpvalue, aux );

      // Skip if we are not calculating this particular value
      if( reduceAtNextStep && !tmpvalue->valueHasBeenSet() ){ 
         members.deactivate(kk); deactivateValue(kk); continue; 
      } else if( !tmpvalue->valueHasBeenSet() ){
         continue;
      }

      // Now incorporate the derivative of the function into the derivatives for the min etc
      keep=false;
      for(unsigned j=0;j<functions.size();++j){
         functions[j]->clear(); functions[j]->calculate( tmpvalue, aux );
         if( functions[j]->sizableContribution( tolerance ) ){ 
             keep=true; functions[j]->mergeDerivatives( kk, *this );
         }
      }  
      // Transfer the value to the field buffers
      if( myfield ){ keep=true; myfield->setBaseQuantity( kk, tmpvalue ); }

      tmpvalue->clearDerivatives();
      for(unsigned i=0;i<aux.size();++i) aux[i].clearDerivatives();
      // If the contribution of this quantity is very small at neighbour list time ignore it
      // untill next neighbour list time
      if( reduceAtNextStep && !keep ){ members.deactivate(kk); deactivateValue(kk); } 
  }
  // Update the dynamic list 
  if(reduceAtNextStep){ members.mpi_gatherActiveMembers( comm ); }
  // MPI Gather everything
  if(!serial){ 
     for(unsigned i=0;i<buffer.size();++i) buffer[i]=0.0;
     unsigned bufsize=0;
     for(unsigned i=0;i<functions.size();++i) functions[i]->copyDataToBuffers( bufsize, buffer ); 
     plumed_assert( bufsize==buffer.size() ); 
     if(buffer.size()>0) comm.Sum( &buffer[0],buffer.size() ); 
     bufsize=0;
     for(unsigned i=0;i<functions.size();++i) functions[i]->retrieveDataFromBuffers( bufsize, buffer ); 
     plumed_assert( bufsize==buffer.size() );

     if( myfield ) myfield->gatherBaseQuantities( comm ); 
  }

  // Set the final value of the function
  for(unsigned j=0;j<final_values.size();++j) functions[j]->finish( final_values[j] ); 

  if( myfield ){
     std::vector<double> thisp( myfield->get_Ndx() ); bool keep;
     Value tmpstress; tmpstress.resizeDerivatives( myfield->get_Ndx() );
     std::vector<Value> tmpder; tmpder.resize( myfield->get_NdX() );
     for(unsigned i=0;i<myfield->get_NdX();++i){
        tmpder[i].set(0);   // This is really important don't delete it
        tmpder[i].resizeDerivatives( myfield->get_Ndx() );
     }
     // Now loop over the spline points
     unsigned ik=0;
     for(unsigned i=0;i<myfield->getNumberOfSplinePoints();++i){
         myfield->getSplinePoint(i,thisp);
         // Calculate the contributions of all the active colvars
         if( updateFreq>0 ){ keep=false; } else { keep=true; }
         for(unsigned j=0;j<members.getNumberActive();++j){
             if( (ik++)%stride!=rank ) continue;  // Ensures we parallelize the double loop over nodes

             kk=members[j];
             //if( myfield->calculateContributionAtPoint( i, kk, this ) ){ keep=true; } 
             unsigned nder=getThisFunctionsNumberOfDerivatives(j);
             if( tmpvalue->getNumberOfDerivatives()!=nder ){ tmpvalue->resizeDerivatives(nder); }
             myfield->extractBaseQuantity( kk, tmpvalue );
             // Calculate the field at point i that arises because of the jth component of the field
             calculateFieldContribution( kk, thisp, tmpvalue, tmpstress, tmpder );
             if( tmpstress.get()>tolerance ){
                 myfield->addStress( i, tmpstress ); keep=true;
                 for(unsigned k=0;k<myfield->get_NdX();++k) myfield->addDerivative( i, k, tmpder[k] );
             }
             // Reset all the tempory values we have used to do this calculation
             tmpvalue->clearDerivatives(); tmpstress.clearDerivatives();
             for(unsigned k=0;k<myfield->get_NdX();++k){ tmpder[k].set(0); tmpder[k].clearDerivatives(); }
         }
         // If the contribution of this quantity is very small at neighbour list time ignore it
         // untill next neighbour list time
         if( reduceAtNextStep && !keep ){ members.deactivate(kk); deactivateValue(kk); }
     }
     // Update the dynamic list 
     if(reduceAtNextStep){ members.mpi_gatherActiveMembers( comm ); }
     // Accumulate the field
     if(!serial) myfield->gatherField( comm );
     // setup the interpolation tables
     myfield->set_tables();
  }
  // Delete the tmpvalues
  delete tmpvalue; 
}

void ActionWithDistribution::retrieveDomain( const unsigned nn, double& min, double& max ){
  plumed_massert(0, "If your function is periodic you need to add a retrieveDomain function so that ActionWithDistribution can retrieve the domain");
}
