/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#include "StoreDataVessel.h"
#include "tools/OpenMP.h"
#include "tools/Stopwatch.h"

using namespace std;
namespace PLMD{
namespace vesselbase{

void ActionWithVessel::registerKeywords(Keywords& keys){
  keys.add("hidden","TOL","this keyword can be used to speed up your calculation. When accumulating sums in which the individual "
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
  keys.addFlag("TIMINGS",false,"output information on the timings of the various parts of the calculation");
  keys.reserveFlag("HIGHMEM",false,"use a more memory intensive version of this collective variable");
  keys.add( vesselRegister().getKeywords() );
}

ActionWithVessel::ActionWithVessel(const ActionOptions&ao):
  Action(ao),
  serial(false),
  lowmem(false),
  noderiv(true),
  actionIsBridged(false),
  mydata(NULL),
  contributorsAreUnlocked(false),
  weightHasDerivatives(false),
  stopwatch(*new Stopwatch)
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
  parseFlag("TIMINGS",timers);
  stopwatch.start(); stopwatch.pause();
}

ActionWithVessel::~ActionWithVessel(){
  for(unsigned i=0;i<functions.size();++i) delete functions[i]; 
  stopwatch.start(); stopwatch.stop();
  if(timers){
     log.printf("timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
     log<<stopwatch;
  }
  delete &stopwatch;
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
  else { delete sv; return; }

  StoreDataVessel* mm=dynamic_cast<StoreDataVessel*>( vv );
  if( mydata && mm ) error("cannot have more than one StoreDataVessel in one action");
  else mydata=mm;
}

BridgeVessel* ActionWithVessel::addBridgingVessel( ActionWithVessel* tome ){
  VesselOptions da("","",0,"",this); 
  BridgeVessel* bv=new BridgeVessel(da);
  bv->setOutputAction( tome );
  tome->actionIsBridged=true;
  functions.push_back( dynamic_cast<Vessel*>(bv) );
  resizeFunctions();
  return bv; 
}

StoreDataVessel* ActionWithVessel::buildDataStashes( const bool& allow_wcutoff, const double& wtol ){
  if(mydata) return mydata;
  
  VesselOptions da("","",0,"",this);
  StoreDataVessel* mm=new StoreDataVessel(da);
  if( allow_wcutoff ) mm->setHardCutoffOnWeight( wtol );
  addVessel(mm);

  // Make sure resizing of vessels is done
  resizeFunctions();

  return mydata;
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
  for(unsigned i=0;i<functions.size();++i) functions[i]->resize();
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

void ActionWithVessel::deactivate_task( const unsigned & task_index ){
  plumed_dbg_assert( contributorsAreUnlocked );
  taskFlags[task_index]=1;
}

void ActionWithVessel::deactivateTasksInRange( const unsigned& lower, const unsigned& upper ){
  plumed_dbg_assert( contributorsAreUnlocked && lower<upper && upper<taskFlags.size() );
  for(unsigned i=lower;i<upper;++i) taskFlags[i]=1;
}

void ActionWithVessel::doJobsRequiredBeforeTaskList(){
  // Do any preparatory stuff for functions
  for(unsigned j=0;j<functions.size();++j) functions[j]->prepare();
}

unsigned ActionWithVessel::getSizeOfBuffer( unsigned& bufsize ){
  for(unsigned i=0;i<functions.size();++i) functions[i]->setBufferStart( bufsize ); 
  if( buffer.size()!=bufsize ) buffer.resize( bufsize );
  if( mydata ){
      unsigned dsize=mydata->getSizeOfDerivativeList();
      if( der_list.size()!=dsize ) der_list.resize( dsize );
  }
  return bufsize;
}

void ActionWithVessel::runAllTasks(){
  if( getExchangeStep() && nactive_tasks!=fullTaskList.size()  ) error("contributors must be unlocked during exchange steps");
  plumed_massert( functions.size()>0, "you must have a call to readVesselKeywords somewhere" );
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){ stride=1; rank=0; }

  // Make sure jobs are done
  if(timers) stopwatch.start("1 Prepare Tasks");
  doJobsRequiredBeforeTaskList();
  if(timers) stopwatch.stop("1 Prepare Tasks");

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nactive_tasks) nt=nactive_tasks/stride/10;
  if( nt==0 ) nt=1;

  // Get size for buffer
  unsigned bsize=0, bufsize=getSizeOfBuffer( bsize ); 
  // Clear buffer
  buffer.assign( buffer.size(), 0.0 );

  // std::vector<unsigned> der_list;
  // if( mydata ) der_list.resize( mydata->getSizeOfDerivativeList(), 0 ); 

  // Build storage stuff for loop
  // std::vector<double> buffer( bufsize, 0.0 );

  if(timers) stopwatch.start("2 Loop over tasks");
#pragma omp parallel num_threads(nt)
{
  std::vector<double> omp_buffer;
  if( nt>1 ) omp_buffer.resize( bufsize, 0.0 );
  MultiValue myvals( getNumberOfQuantities(), getNumberOfDerivatives() );
  MultiValue bvals( getNumberOfQuantities(), getNumberOfDerivatives() );
  myvals.clearAll(); bvals.clearAll();
 
#pragma omp for nowait
  for(unsigned i=rank;i<nactive_tasks;i+=stride){
      // Calculate the stuff in the loop for this action
      performTask( indexOfTaskInFullList[i], partialTaskList[i], myvals );
      // Weight should be between zero and one
      plumed_dbg_assert( myvals.get(0)>=0 && myvals.get(0)<=1.0 );

      // Check for conditions that allow us to just to skip the calculation
      // the condition is that the weight of the contribution is low 
      // N.B. Here weights are assumed to be between zero and one
      if( myvals.get(0)<tolerance ){
         // Deactivate task if it is less than the neighbor list tolerance
         if( myvals.get(0)<nl_tolerance && contributorsAreUnlocked ) deactivate_task( indexOfTaskInFullList[i] );
         // Clear the derivatives
         myvals.clearAll();
         continue;
      }

      // Now calculate all the functions
      // If the contribution of this quantity is very small at neighbour list time ignore it
      // untill next neighbour list time
      if( nt>1 ){
          if( !calculateAllVessels( indexOfTaskInFullList[i], myvals, bvals, omp_buffer, der_list ) && contributorsAreUnlocked ) deactivate_task( indexOfTaskInFullList[i] );
      } else {
          if( !calculateAllVessels( indexOfTaskInFullList[i], myvals, bvals, buffer, der_list ) && contributorsAreUnlocked ) deactivate_task( indexOfTaskInFullList[i] );
      }

      // Clear the value
      myvals.clearAll();
  }
#pragma omp critical
  if(nt>1) for(unsigned i=0;i<bufsize;++i) buffer[i]+=omp_buffer[i];
}
  if(timers) stopwatch.stop("2 Loop over tasks");

  if(timers) stopwatch.start("3 MPI gather");
  // MPI Gather everything
  if( !serial && buffer.size()>0 ) comm.Sum( buffer );
  // MPI Gather index stores
  if( mydata && !lowmem && !noderiv ){ 
     comm.Sum( der_list ); mydata->setActiveValsAndDerivatives( der_list ); 
  }
  // Update the elements that are makign contributions to the sum here
  // this causes problems if we do it in prepare
  if( !serial && contributorsAreUnlocked ) comm.Sum( taskFlags );
  if(timers) stopwatch.stop("3 MPI gather");

  if(timers) stopwatch.start("4 Finishing computations");
  finishComputations( buffer );
  if(timers) stopwatch.stop("4 Finishing computations");
}

void ActionWithVessel::transformBridgedDerivatives( const unsigned& current, MultiValue& invals, MultiValue& outvals ) const {
  plumed_error();
}

bool ActionWithVessel::calculateAllVessels( const unsigned& taskCode, MultiValue& myvals, MultiValue& bvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ){
  bool keep=false; 
  for(unsigned j=0;j<functions.size();++j){
      // Calculate returns a bool that tells us if this particular
      // quantity is contributing more than the tolerance
      if( functions[j]->calculate( taskCode, functions[j]->transformDerivatives(taskCode, myvals, bvals), buffer, der_list ) ) keep=true;
      if( !actionIsBridged && bvals.getNumberActive()>0 ) bvals.clearAll(); 
  }
  return keep;
}

void ActionWithVessel::finishComputations( const std::vector<double>& buffer ){
  // Set the final value of the function
  for(unsigned j=0;j<functions.size();++j) functions[j]->finish( buffer );  
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
