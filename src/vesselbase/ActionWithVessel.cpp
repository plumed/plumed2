/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#include "tools/Communicator.h"
#include "ActionWithVessel.h"
#include "Vessel.h"
#include "ShortcutVessel.h"
#include "StoreDataVessel.h"
#include "VesselRegister.h"
#include "BridgeVessel.h"
#include "FunctionVessel.h"

using namespace std;
namespace PLMD{
namespace vesselbase{

void ActionWithVessel::registerKeywords(Keywords& keys){
  keys.add("optional","TOL","this keyword can be used to speed up your calculation. When accumulating sums in which the individual "
                            "terms are numbers inbetween zero and one it is assumed that terms less than a certain tolerance "
                            "make only a small contribution to the sum.  They can thus be safely ignored as can the the derivatives "
                            "wrt these small quantities.");
  keys.reserve("hidden","NL_TOL","this keyword can be used to speed up your calculation.  It must be used in conjuction with the TOL "
                                   "keyword and the value for NL_TOL must be set less than the value for TOL.  This keyword ensures that "
                                   "quantities, which are much less than TOL and which will thus not added to the sums being accumulated "
                                   "are not calculated at every step. They are only calculated when the neighbor list is updated.");
  keys.add("hidden","MAXDERIVATIVES","The maximum number of derivatives that can be used when storing data.  This controls when "
                                     "we have to start using lowmem");
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize");
  keys.addFlag("LOWMEM",false,"lower the memory requirements");
  keys.reserveFlag("HIGHMEM",false,"use a more memory intensive version of this collective variable");
  keys.add( vesselRegister().getKeywords() );
}

ActionWithVessel::ActionWithVessel(const ActionOptions&ao):
  Action(ao),
  serial(false),
  lowmem(false),
  noderiv(true),
  contributorsAreUnlocked(false),
  weightHasDerivatives(false)
{
  maxderivatives=309; parse("MAXDERIVATIVES",maxderivatives);
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  else serial=true;
  if(serial)log.printf("  doing calculation in serial\n");
  if( keywords.exists("LOWMEM") ){
     plumed_assert( !keywords.exists("HIGHMEM") );
     parseFlag("LOWMEM",lowmem);
     if(lowmem)log.printf("  lowering memory requirements\n");
  } 
  if( keywords.exists("HIGHMEM") ){
     plumed_assert( !keywords.exists("LOWMEM") );
     bool highmem; parseFlag("HIGHMEM",highmem);
     lowmem=!highmem;
     if(!lowmem) log.printf("  increasing the memory requirements\n");
  }
  tolerance=nl_tolerance=epsilon; 
  if( keywords.exists("TOL") ) parse("TOL",tolerance);
  if( tolerance>epsilon){
     if( keywords.exists("NL_TOL") ) parse("NL_TOL",nl_tolerance);
     if( nl_tolerance>tolerance ) error("NL_TOL must be smaller than TOL"); 
     log.printf(" Ignoring contributions less than %f",tolerance);
     if( nl_tolerance>epsilon ) log.printf(" and ignoring quantities less than %f inbetween neighbor list update steps\n",nl_tolerance);
     else log.printf("\n");
  }
}

ActionWithVessel::~ActionWithVessel(){
  for(unsigned i=0;i<functions.size();++i) delete functions[i]; 
}

void ActionWithVessel::addVessel( const std::string& name, const std::string& input, const int numlab ){
  VesselOptions da(name,"",numlab,input,this);
  Vessel* vv=vesselRegister().create(name,da);
  FunctionVessel* fv=dynamic_cast<FunctionVessel*>(vv);
  if( fv ){
      std::string mylabel=Vessel::transformName( name );
      plumed_massert( keywords.outputComponentExists(mylabel,false), "a description of the value calculated by vessel " + name + " has not been added to the manual"); 
  } 
  addVessel(vv);
}

void ActionWithVessel::addVessel( Vessel* vv ){
  ShortcutVessel* sv=dynamic_cast<ShortcutVessel*>(vv);
  if(!sv){ vv->checkRead(); functions.push_back(vv); }
  else { delete sv; }
}

BridgeVessel* ActionWithVessel::addBridgingVessel( ActionWithVessel* tome ){
  VesselOptions da("","",0,"",this); 
  BridgeVessel* bv=new BridgeVessel(da);
  bv->setOutputAction( tome );
  functions.push_back( dynamic_cast<Vessel*>(bv) );
  resizeFunctions();
  return bv; 
}

StoreDataVessel* ActionWithVessel::buildDataStashes( const bool& allow_wcutoff, const double& wtol ){
  for(unsigned i=0;i<functions.size();++i){
      StoreDataVessel* vsv=dynamic_cast<StoreDataVessel*>( functions[i] );
      if( vsv ) return vsv;
  }
  return NULL;
}

void ActionWithVessel::addTaskToList( const unsigned& taskCode ){
  indexOfTaskInFullList.push_back( fullTaskList.size() );
  fullTaskList.push_back( taskCode ); partialTaskList.push_back( taskCode ); 
  taskFlags.push_back(0); nactive_tasks = fullTaskList.size();
  plumed_assert( partialTaskList.size()==nactive_tasks && indexOfTaskInFullList.size()==nactive_tasks && taskFlags.size()==nactive_tasks );
}

void ActionWithVessel::readVesselKeywords(){
  // Set maxderivatives if it is too big
  if( maxderivatives>getNumberOfDerivatives() ) maxderivatives=getNumberOfDerivatives();

  // Loop over all keywords find the vessels and create appropriate functions
  for(unsigned i=0;i<keywords.size();++i){
      std::string thiskey,input; thiskey=keywords.getKeyword(i);
      // Check if this is a key for a vessel
      if( vesselRegister().check(thiskey) ){
          // If the keyword is a flag read it in as a flag
          if( keywords.style(thiskey,"flag") ){
              bool dothis; parseFlag(thiskey,dothis);
              if(dothis) addVessel( thiskey, input );
          // If it is numbered read it in as a numbered thing
          } else if( keywords.numbered(thiskey) ) {
              parse(thiskey,input);
              if(input.size()!=0){ 
                    addVessel( thiskey, input );
              } else {
                 for(unsigned i=1;;++i){
                    if( !parseNumbered(thiskey,i,input) ) break;
                    std::string ss; Tools::convert(i,ss);
                    addVessel( thiskey, input, i ); 
                    input.clear();
                 } 
              }
          // Otherwise read in the keyword the normal way
          } else {
              parse(thiskey, input);
              if(input.size()!=0) addVessel(thiskey,input);
          }
          input.clear();
      }
  }

  // Make sure all vessels have had been resized at start
  if( functions.size()>0 ) resizeFunctions();
}

void ActionWithVessel::resizeFunctions(){
  unsigned bufsize=0; 
  for(unsigned i=0;i<functions.size();++i){
     functions[i]->bufstart=bufsize;
     functions[i]->resize();
     bufsize+=functions[i]->bufsize;
  }
  thisval.resize( getNumberOfQuantities() ); thisval_wasset.resize( getNumberOfQuantities(), false );
  derivatives.resize( getNumberOfQuantities()*getNumberOfDerivatives(), 0.0 );
  buffer.resize( bufsize );
}

void ActionWithVessel::needsDerivatives(){
  // Turn on the derivatives and resize
  noderiv=false; resizeFunctions(); 
  // Setting contributors unlocked here ensures that link cells are ignored
  contributorsAreUnlocked=true; finishTaskListUpdate(); contributorsAreUnlocked=false;
  // And turn on the derivatives in all actions on which we are dependent
  for(unsigned i=0;i<getDependencies().size();++i){
      ActionWithVessel* vv=dynamic_cast<ActionWithVessel*>( getDependencies()[i] );
      if(vv) vv->needsDerivatives();
  }
}

void ActionWithVessel::unlockContributors(){
  if( contributorsAreUnlocked ) return;
  nactive_tasks = fullTaskList.size();
  for(unsigned i=0;i<fullTaskList.size();++i){ 
     partialTaskList[i] = fullTaskList[i]; taskFlags[i]=0; 
     indexOfTaskInFullList[i]=i;
  }
  finishTaskListUpdate();
  contributorsAreUnlocked=true;
  resizeFunctions();
}

void ActionWithVessel::lockContributors(){
  nactive_tasks = 0;
  for(unsigned i=0;i<fullTaskList.size();++i){
      // Deactivate sets inactive tasks to number not equal to zero
      if( taskFlags[i]==0 ){
          partialTaskList[nactive_tasks] = fullTaskList[i]; 
          indexOfTaskInFullList[nactive_tasks]=i;
          nactive_tasks++; 
      } 
  }
  contributorsAreUnlocked=false;
  finishTaskListUpdate();
  resizeFunctions();
}

void ActionWithVessel::deactivateAllTasks(){
  contributorsAreUnlocked=true; nactive_tasks = 0;
}

void ActionWithVessel::activateTheseTasks( std::vector<unsigned>& additionalTasks ){
  plumed_dbg_assert( additionalTasks.size()==fullTaskList.size() );
  // Activate tasks that are already active locally
  for(unsigned i=0;i<nactive_tasks;++i) additionalTasks[ indexOfTaskInFullList[i] ] = 1;

  nactive_tasks = 0;
  for(unsigned i=0;i<fullTaskList.size();++i){
      // Deactivate sets inactive tasks to number not equal to zero
      if( additionalTasks[i]>0 ){
          partialTaskList[nactive_tasks] = fullTaskList[i]; 
          indexOfTaskInFullList[nactive_tasks]=i;
          nactive_tasks++;
      } else {
          taskFlags[i]=1;
      }
  }
  contributorsAreUnlocked=false;
}

void ActionWithVessel::deactivate_task(){
  plumed_dbg_assert( contributorsAreUnlocked );
  taskFlags[task_index]=1;
}

void ActionWithVessel::deactivateTasksInRange( const unsigned& lower, const unsigned& upper ){
  plumed_dbg_assert( contributorsAreUnlocked && lower<upper && upper<taskFlags.size() );
  for(unsigned i=lower;i<upper;++i) taskFlags[i]=1;
}

void ActionWithVessel::doJobsRequiredBeforeTaskList(){
  // Clear all data from previous calculations
  buffer.assign(buffer.size(),0.0);
  // Do any preparatory stuff for functions
  for(unsigned j=0;j<functions.size();++j) functions[j]->prepare();
}

void ActionWithVessel::runAllTasks(){
  if( getExchangeStep() && nactive_tasks!=fullTaskList.size()  ) error("contributors must be unlocked during exchange steps");
  plumed_massert( functions.size()>0, "you must have a call to readVesselKeywords somewhere" );
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){ stride=1; rank=0; }

  // Make sure jobs are done
  doJobsRequiredBeforeTaskList();

  for(unsigned i=rank;i<nactive_tasks;i+=stride){
      // The index of the task in the full list
      task_index=indexOfTaskInFullList[i];
      // Store the task we are currently working on
      current=partialTaskList[i];
      // Calculate the stuff in the loop for this action
      performTask();
      // Weight should be between zero and one
      plumed_dbg_assert( getValueForTolerance()>=0 && getValueForTolerance()<=1.0 );

      // Check for conditions that allow us to just to skip the calculation
      // the condition is that the weight of the contribution is low 
      // N.B. Here weights are assumed to be between zero and one
      if( getValueForTolerance()<tolerance ){
         // Clear the derivatives
         clearAfterTask();  
         // Deactivate task if it is less than the neighbor list tolerance
         if( getValueForTolerance()<nl_tolerance && contributorsAreUnlocked ) deactivate_task();
         continue;
      }

      // Now calculate all the functions
      // If the contribution of this quantity is very small at neighbour list time ignore it
      // untill next neighbour list time
      if( !calculateAllVessels() && contributorsAreUnlocked ) deactivate_task();
  }
  finishComputations();
}

void ActionWithVessel::getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices ){
  indices[jstore]=getNumberOfDerivatives();
  if( indices[jstore]>maxder ) error("too many derivatives to store. Run with LOWMEM");

  unsigned kder = ntotal + jstore*getNumberOfDerivatives();
  for(unsigned jder=0;jder<getNumberOfDerivatives();++jder){ indices[ kder ] = jder; kder++; }
}

void ActionWithVessel::clearAfterTask(){
  // Clear the derivatives from this step
  for(unsigned k=0;k<thisval.size();++k) clearDerivativesAfterTask(k);
}

void ActionWithVessel::clearDerivativesAfterTask( const unsigned& ider ){
  thisval[ider]=0.0; thisval_wasset[ider]=false;
  if( !noderiv ){
     unsigned kstart=ider*getNumberOfDerivatives();
     for(unsigned j=0;j<getNumberOfDerivatives();++j) derivatives[ kstart+j ]=0.0;
  }
}

bool ActionWithVessel::calculateAllVessels(){
  bool keep=false;
  for(unsigned j=0;j<functions.size();++j){
      // Calculate returns a bool that tells us if this particular
      // quantity is contributing more than the tolerance
      if( functions[j]->calculate() ) keep=true;
  }
  clearAfterTask();
  return keep;
}

void ActionWithVessel::finishComputations(){
  // MPI Gather everything
  if( !serial && buffer.size()>0 ) comm.Sum( buffer );
  // Update the elements that are makign contributions to the sum here
  // this causes problems if we do it in prepare
  if( !serial && contributorsAreUnlocked ) comm.Sum( taskFlags );

  // Set the final value of the function
  for(unsigned j=0;j<functions.size();++j) functions[j]->finish(); 
}

void ActionWithVessel::chainRuleForElementDerivatives( const unsigned& iout, const unsigned& ider, const double& df, Vessel* valout ){
  if( noderiv ) return;
  current_buffer_stride=1;
  current_buffer_start=valout->bufstart + (getNumberOfDerivatives()+1)*iout + 1;
  mergeDerivatives( ider, df );
} 

void ActionWithVessel::chainRuleForElementDerivatives( const unsigned& iout, const unsigned& ider, const unsigned& stride, 
                                                       const unsigned& off, const double& df, Vessel* valout ){
  if( noderiv ) return;
  plumed_dbg_assert( off<stride );
  current_buffer_stride=stride;
  current_buffer_start=valout->bufstart + stride*(getNumberOfDerivatives()+1)*iout + stride + off;
  mergeDerivatives( ider, df );
}

void ActionWithVessel::mergeDerivatives( const unsigned& ider, const double& df ){
  unsigned nder=getNumberOfDerivatives(), vstart=nder*ider; 
  for(unsigned i=0;i<getNumberOfDerivatives();++i){
     accumulateDerivative( i, df*derivatives[vstart+i] ); 
  }
}

bool ActionWithVessel::getForcesFromVessels( std::vector<double>& forcesToApply ){
#ifndef DNDEBUG
  if( forcesToApply.size()>0 ) plumed_dbg_assert( forcesToApply.size()==getNumberOfDerivatives() );
#endif
  if(tmpforces.size()!=forcesToApply.size() ) tmpforces.resize( forcesToApply.size() );

  forcesToApply.assign( forcesToApply.size(),0.0 );
  bool wasforced=false;
  for(unsigned i=0;i<getNumberOfVessels();++i){
    if( (functions[i]->applyForce( tmpforces )) ){
       wasforced=true;
       for(unsigned j=0;j<forcesToApply.size();++j) forcesToApply[j]+=tmpforces[j];
    }
  }
  return wasforced;
}

void ActionWithVessel::retrieveDomain( std::string& min, std::string& max ){
  plumed_merror("If your function is periodic you need to add a retrieveDomain function so that ActionWithVessel can retrieve the domain");
}

Vessel* ActionWithVessel::getVesselWithName( const std::string& mynam ){
  int target=-1;
  for(unsigned i=0;i<functions.size();++i){
     if( functions[i]->getName().find(mynam)!=std::string::npos ){
        if( target<0 ) target=i; 
        else error("found more than one " + mynam + " object in action");
     }  
  }
  return functions[target];
}

}
}
