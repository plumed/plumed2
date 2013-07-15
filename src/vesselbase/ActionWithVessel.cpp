/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "tools/Communicator.h"
#include "ActionWithVessel.h"
#include "Vessel.h"
#include "ShortcutVessel.h"
#include "VesselRegister.h"
#include "BridgeVessel.h"

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
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize");
  keys.addFlag("LOWMEM",false,"lower the memory requirements");
  keys.add( vesselRegister().getKeywords() );
}

ActionWithVessel::ActionWithVessel(const ActionOptions&ao):
  Action(ao),
  read(false),
  serial(false),
  contributorsAreUnlocked(false),
  weightHasDerivatives(false)
{
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  else serial=true;
  if(serial)log.printf("  doing calculation in serial\n");
  parseFlag("LOWMEM",lowmem);
  if(lowmem)log.printf("  lowering memory requirements\n");
  tolerance=nl_tolerance=epsilon; 
  if( keywords.exists("TOL") ) parse("TOL",tolerance);
  if( tolerance>epsilon){
     if( keywords.exists("NL_TOL") ) parse("NL_TOL",nl_tolerance);
     if( nl_tolerance>tolerance ) error("NL_TOL must be smaller than TOL"); 
     log.printf(" Ignoring contributions less than %lf",tolerance);
     if( nl_tolerance>epsilon ) log.printf(" and ignoring quantities less than %lf inbetween neighbor list update steps\n",nl_tolerance);
     else log.printf("\n");
  }
  // Setup stuff for communicating what tasks have been deactivated across all nodes
  taskList.setupMPICommunication( comm );
}

ActionWithVessel::~ActionWithVessel(){
  for(unsigned i=0;i<functions.size();++i) delete functions[i]; 
}

void ActionWithVessel::addVessel( const std::string& name, const std::string& input, const int numlab, const std::string thislab ){
  read=true; VesselOptions da(name,thislab,numlab,input,this);
  Vessel* vv=vesselRegister().create(name,da); vv->checkRead();
  addVessel(vv);
}

void ActionWithVessel::addVessel( Vessel* vv ){
  ShortcutVessel* sv=dynamic_cast<ShortcutVessel*>(vv);
  if(!sv) functions.push_back(vv);
}

BridgeVessel* ActionWithVessel::addBridgingVessel( ActionWithVessel* tome ){
  read=true; VesselOptions da("","",0,"",this); 
  BridgeVessel* bv=new BridgeVessel(da);
  bv->setOutputAction( tome );
  functions.push_back( dynamic_cast<Vessel*>(bv) );
  resizeFunctions();
  return bv; 
}

void ActionWithVessel::readVesselKeywords(){
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
  tmpforces.resize( getNumberOfDerivatives() );
}

void ActionWithVessel::deactivate_task(){
  plumed_dbg_assert( contributorsAreUnlocked );
  for(unsigned i=0;i<taskList.fullSize();++i){
     if( taskList(i)==current ) taskList.deactivate(i);
  }
}

void ActionWithVessel::activateTheseTasks( const std::vector<bool>& additionalTasks ){
  plumed_dbg_assert( additionalTasks.size()==taskList.fullSize() );
  for(unsigned i=0;i<additionalTasks.size();++i){
      if( additionalTasks[i] ) taskList.activate(i);
  }
  taskList.updateActiveMembers();
}

void ActionWithVessel::doJobsRequiredBeforeTaskList(){
  // Clear all data from previous calculations
  buffer.assign(buffer.size(),0.0);
  // Do any preparatory stuff for functions
  for(unsigned j=0;j<functions.size();++j) functions[j]->prepare();
}

void ActionWithVessel::runAllTasks(){
  plumed_massert( read, "you must have a call to readVesselKeywords somewhere" );
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){ stride=1; rank=0; }

  // Make sure jobs are done
  doJobsRequiredBeforeTaskList();

  for(unsigned i=rank;i<taskList.getNumberActive();i+=stride){
      // Store the task we are currently working on
      current=taskList[i];
      // Calculate the stuff in the loop for this action
      performTask(i);
      // Weight should be between zero and one
      plumed_dbg_assert( thisval[1]>=0 && thisval[1]<=1.0 );

      // Check for conditions that allow us to just to skip the calculation
      // the condition is that the weight of the contribution is low 
      // N.B. Here weights are assumed to be between zero and one
      if( thisval[1]<tolerance ){
         // Clear the derivatives
         clearAfterTask();  
         // Deactivate task if it is less than the neighbor list tolerance
         if( thisval[1]<nl_tolerance && contributorsAreUnlocked ) deactivate_task();
         continue;
      }

      // Now calculate all the functions
      // If the contribution of this quantity is very small at neighbour list time ignore it
      // untill next neighbour list time
      if( !calculateAllVessels() && contributorsAreUnlocked ) deactivate_task();
  }
  finishComputations();
}

void ActionWithVessel::clearAfterTask(){
  // Clear the derivatives from this step
  for(unsigned k=0;k<thisval.size();++k){
     thisval_wasset[k]=false;
     clearDerivativesAfterTask(k);
  }
}

void ActionWithVessel::clearDerivativesAfterTask( const unsigned& ider ){
  unsigned kstart=ider*getNumberOfDerivatives();
  for(unsigned j=0;j<getNumberOfDerivatives();++j) derivatives[ kstart+j ]=0.0;
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
  if(!serial && buffer.size()>0) comm.Sum( buffer);

  // Set the final value of the function
  for(unsigned j=0;j<functions.size();++j) functions[j]->finish(); 
}

void ActionWithVessel::chainRuleForElementDerivatives( const unsigned& iout, const unsigned& ider, const double& df, Vessel* valout ){
  current_buffer_stride=1;
  current_buffer_start=valout->bufstart + (getNumberOfDerivatives()+1)*iout + 1;
  mergeDerivatives( ider, df );
} 

void ActionWithVessel::chainRuleForElementDerivatives( const unsigned& iout, const unsigned& ider, const unsigned& stride, 
                                                       const unsigned& off, const double& df, Vessel* valout ){
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
  plumed_dbg_assert( forcesToApply.size()==getNumberOfDerivatives() );
  forcesToApply.assign( forcesToApply.size(),0.0 );
  bool wasforced=false;
  for(int i=0;i<getNumberOfVessels();++i){
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

}
}
