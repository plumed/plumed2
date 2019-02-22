/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "ActionWithVessel.h"
#include "tools/Communicator.h"
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
namespace PLMD {
namespace vesselbase {

void ActionWithVessel::registerKeywords(Keywords& keys) {
  keys.add("hidden","TOL","this keyword can be used to speed up your calculation. When accumulating sums in which the individual "
           "terms are numbers in between zero and one it is assumed that terms less than a certain tolerance "
           "make only a small contribution to the sum.  They can thus be safely ignored as can the the derivatives "
           "wrt these small quantities.");
  keys.add("hidden","MAXDERIVATIVES","The maximum number of derivatives that can be used when storing data.  This controls when "
           "we have to start using lowmem");
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not use MPI");
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
  nactive_tasks(0),
  dertime_can_be_off(false),
  dertime(true),
  contributorsAreUnlocked(false),
  weightHasDerivatives(false),
  mydata(NULL)
{
  maxderivatives=309; parse("MAXDERIVATIVES",maxderivatives);
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  else serial=true;
  if(serial)log.printf("  doing calculation in serial\n");
  if( keywords.exists("LOWMEM") ) {
    plumed_assert( !keywords.exists("HIGHMEM") );
    parseFlag("LOWMEM",lowmem);
    if(lowmem) {
      log.printf("  lowering memory requirements\n");
      dertime_can_be_off=true;
    }
  }
  if( keywords.exists("HIGHMEM") ) {
    plumed_assert( !keywords.exists("LOWMEM") );
    bool highmem; parseFlag("HIGHMEM",highmem);
    lowmem=!highmem;
    if(!lowmem) log.printf("  increasing the memory requirements\n");
  }
  tolerance=nl_tolerance=epsilon;
  if( keywords.exists("TOL") ) parse("TOL",tolerance);
  if( tolerance>epsilon) {
    log.printf(" Ignoring contributions less than %f \n",tolerance);
  }
  parseFlag("TIMINGS",timers);
  stopwatch.start(); stopwatch.pause();
}

ActionWithVessel::~ActionWithVessel() {
  stopwatch.start(); stopwatch.stop();
  if(timers) {
    log.printf("timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log<<stopwatch;
  }
}

void ActionWithVessel::addVessel( const std::string& name, const std::string& input, const int numlab ) {
  VesselOptions da(name,"",numlab,input,this);
  auto vv=vesselRegister().create(name,da);
  FunctionVessel* fv=dynamic_cast<FunctionVessel*>(vv.get());
  if( fv ) {
    std::string mylabel=Vessel::transformName( name );
    plumed_massert( keywords.outputComponentExists(mylabel,false), "a description of the value calculated by vessel " + name + " has not been added to the manual");
  }
  addVessel(std::move(vv));
}

void ActionWithVessel::addVessel( std::unique_ptr<Vessel> vv_ptr ) {

// In the original code, the dynamically casted pointer was deleted here.
// Now that vv_ptr is a unique_ptr, the object will be deleted automatically when
// exiting this routine.
  if(dynamic_cast<ShortcutVessel*>(vv_ptr.get())) return;

  vv_ptr->checkRead();

  StoreDataVessel* mm=dynamic_cast<StoreDataVessel*>( vv_ptr.get() );
  if( mydata && mm ) error("cannot have more than one StoreDataVessel in one action");
  else if( mm ) mydata=mm;
  else dertime_can_be_off=false;

// Ownership is transfered to functions
  functions.emplace_back(std::move(vv_ptr));
}

BridgeVessel* ActionWithVessel::addBridgingVessel( ActionWithVessel* tome ) {
  VesselOptions da("","",0,"",this);
  std::unique_ptr<BridgeVessel> bv(new BridgeVessel(da));
  bv->setOutputAction( tome );
  tome->actionIsBridged=true; dertime_can_be_off=false;
// store this pointer in order to return it later.
// notice that I cannot access this with functions.tail().get()
// since functions contains pointers to a different class (Vessel)
  auto toBeReturned=bv.get();
  functions.emplace_back( std::move(bv) );
  resizeFunctions();
  return toBeReturned;
}

StoreDataVessel* ActionWithVessel::buildDataStashes( ActionWithVessel* actionThatUses ) {
  if(mydata) {
    if( actionThatUses ) mydata->addActionThatUses( actionThatUses );
    return mydata;
  }

  VesselOptions da("","",0,"",this);
  std::unique_ptr<StoreDataVessel> mm( new StoreDataVessel(da) );
  if( actionThatUses ) mm->addActionThatUses( actionThatUses );
  addVessel(std::move(mm));

  // Make sure resizing of vessels is done
  resizeFunctions();

  return mydata;
}

void ActionWithVessel::addTaskToList( const unsigned& taskCode ) {
  fullTaskList.push_back( taskCode ); taskFlags.push_back(0);
  plumed_assert( fullTaskList.size()==taskFlags.size() );
}

void ActionWithVessel::readVesselKeywords() {
  // Set maxderivatives if it is too big
  if( maxderivatives>getNumberOfDerivatives() ) maxderivatives=getNumberOfDerivatives();

  // Loop over all keywords find the vessels and create appropriate functions
  for(unsigned i=0; i<keywords.size(); ++i) {
    std::string thiskey,input; thiskey=keywords.getKeyword(i);
    // Check if this is a key for a vessel
    if( vesselRegister().check(thiskey) ) {
      plumed_assert( keywords.style(thiskey,"vessel") );
      bool dothis=false; parseFlag(thiskey,dothis);
      if(dothis) addVessel( thiskey, input );

      parse(thiskey,input);
      if(input.size()!=0) {
        addVessel( thiskey, input );
      } else {
        for(unsigned i=1;; ++i) {
          if( !parseNumbered(thiskey,i,input) ) break;
          std::string ss; Tools::convert(i,ss);
          addVessel( thiskey, input, i );
          input.clear();
        }
      }
    }
  }

  // Make sure all vessels have had been resized at start
  if( functions.size()>0 ) resizeFunctions();
}

void ActionWithVessel::resizeFunctions() {
  for(unsigned i=0; i<functions.size(); ++i) functions[i]->resize();
}

void ActionWithVessel::needsDerivatives() {
  // Turn on the derivatives and resize
  noderiv=false; resizeFunctions();
  // Setting contributors unlocked here ensures that link cells are ignored
  contributorsAreUnlocked=true; contributorsAreUnlocked=false;
  // And turn on the derivatives in all actions on which we are dependent
  for(unsigned i=0; i<getDependencies().size(); ++i) {
    ActionWithVessel* vv=dynamic_cast<ActionWithVessel*>( getDependencies()[i] );
    if(vv) vv->needsDerivatives();
  }
}

void ActionWithVessel::lockContributors() {
  nactive_tasks = 0;
  for(unsigned i=0; i<fullTaskList.size(); ++i) {
    if( taskFlags[i]>0 ) nactive_tasks++;
  }

  unsigned n=0;
  partialTaskList.resize( nactive_tasks );
  indexOfTaskInFullList.resize( nactive_tasks );
  for(unsigned i=0; i<fullTaskList.size(); ++i) {
    // Deactivate sets inactive tasks to number not equal to zero
    if( taskFlags[i]>0 ) {
      partialTaskList[n] = fullTaskList[i];
      indexOfTaskInFullList[n]=i;
      n++;
    }
  }
  plumed_dbg_assert( n==nactive_tasks );
  for(unsigned i=0; i<functions.size(); ++i) {
    BridgeVessel* bb = dynamic_cast<BridgeVessel*>( functions[i].get() );
    if( bb ) bb->copyTaskFlags();
  }
  // Resize mydata to accomodate all active tasks
  if( mydata ) mydata->resize();
  contributorsAreUnlocked=false;
}

void ActionWithVessel::deactivateAllTasks() {
  contributorsAreUnlocked=true; nactive_tasks = 0;
  taskFlags.assign(taskFlags.size(),0);
}

bool ActionWithVessel::taskIsCurrentlyActive( const unsigned& index ) const {
  plumed_dbg_assert( index<taskFlags.size() ); return (taskFlags[index]>0);
}

void ActionWithVessel::doJobsRequiredBeforeTaskList() {
  // Do any preparatory stuff for functions
  for(unsigned j=0; j<functions.size(); ++j) functions[j]->prepare();
}

unsigned ActionWithVessel::getSizeOfBuffer( unsigned& bufsize ) {
  for(unsigned i=0; i<functions.size(); ++i) functions[i]->setBufferStart( bufsize );
  if( buffer.size()!=bufsize ) buffer.resize( bufsize );
  if( mydata ) {
    unsigned dsize=mydata->getSizeOfDerivativeList();
    if( der_list.size()!=dsize ) der_list.resize( dsize );
  }
  return bufsize;
}

void ActionWithVessel::runAllTasks() {
  plumed_massert( !contributorsAreUnlocked && functions.size()>0, "you must have a call to readVesselKeywords somewhere" );
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) { stride=1; rank=0; }

  // Make sure jobs are done
  if(timers) stopwatch.start("1 Prepare Tasks");
  doJobsRequiredBeforeTaskList();
  if(timers) stopwatch.stop("1 Prepare Tasks");

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*2>nactive_tasks || !threadSafe()) nt=1;

  // Get size for buffer
  unsigned bsize=0, bufsize=getSizeOfBuffer( bsize );
  // Clear buffer
  buffer.assign( buffer.size(), 0.0 );
  // Switch off calculation of derivatives in main loop
  if( dertime_can_be_off ) dertime=false;

  if(timers) stopwatch.start("2 Loop over tasks");
  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_buffer;
    if( nt>1 ) omp_buffer.resize( bufsize, 0.0 );
    MultiValue myvals( getNumberOfQuantities(), getNumberOfDerivatives() );
    MultiValue bvals( getNumberOfQuantities(), getNumberOfDerivatives() );
    myvals.clearAll(); bvals.clearAll();

    #pragma omp for nowait schedule(dynamic)
    for(unsigned i=rank; i<nactive_tasks; i+=stride) {
      // Calculate the stuff in the loop for this action
      performTask( indexOfTaskInFullList[i], partialTaskList[i], myvals );

      // Check for conditions that allow us to just to skip the calculation
      // the condition is that the weight of the contribution is low
      // N.B. Here weights are assumed to be between zero and one
      if( myvals.get(0)<tolerance ) {
        // Clear the derivatives
        myvals.clearAll();
        continue;
      }

      // Now calculate all the functions
      // If the contribution of this quantity is very small at neighbour list time ignore it
      // untill next neighbour list time
      if( nt>1 ) {
        calculateAllVessels( indexOfTaskInFullList[i], myvals, bvals, omp_buffer, der_list );
      } else {
        calculateAllVessels( indexOfTaskInFullList[i], myvals, bvals, buffer, der_list );
      }

      // Clear the value
      myvals.clearAll();
    }
    #pragma omp critical
    if(nt>1) for(unsigned i=0; i<bufsize; ++i) buffer[i]+=omp_buffer[i];
  }
  if(timers) stopwatch.stop("2 Loop over tasks");
  // Turn back on derivative calculation
  dertime=true;

  if(timers) stopwatch.start("3 MPI gather");
  // MPI Gather everything
  if( !serial && buffer.size()>0 ) comm.Sum( buffer );
  // MPI Gather index stores
  if( mydata && !lowmem && !noderiv ) {
    comm.Sum( der_list ); mydata->setActiveValsAndDerivatives( der_list );
  }
  // Update the elements that are makign contributions to the sum here
  // this causes problems if we do it in prepare
  if(timers) stopwatch.stop("3 MPI gather");

  if(timers) stopwatch.start("4 Finishing computations");
  finishComputations( buffer );
  if(timers) stopwatch.stop("4 Finishing computations");
}

void ActionWithVessel::transformBridgedDerivatives( const unsigned& current, MultiValue& invals, MultiValue& outvals ) const {
  plumed_error();
}

void ActionWithVessel::calculateAllVessels( const unsigned& taskCode, MultiValue& myvals, MultiValue& bvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) {
  for(unsigned j=0; j<functions.size(); ++j) {
    // Calculate returns a bool that tells us if this particular
    // quantity is contributing more than the tolerance
    functions[j]->calculate( taskCode, functions[j]->transformDerivatives(taskCode, myvals, bvals), buffer, der_list );
    if( !actionIsBridged ) bvals.clearAll();
  }
  return;
}

void ActionWithVessel::finishComputations( const std::vector<double>& buffer ) {
  // Set the final value of the function
  for(unsigned j=0; j<functions.size(); ++j) functions[j]->finish( buffer );
}

bool ActionWithVessel::getForcesFromVessels( std::vector<double>& forcesToApply ) {
#ifndef NDEBUG
  if( forcesToApply.size()>0 ) plumed_dbg_assert( forcesToApply.size()==getNumberOfDerivatives() );
#endif
  if(tmpforces.size()!=forcesToApply.size() ) tmpforces.resize( forcesToApply.size() );

  forcesToApply.assign( forcesToApply.size(),0.0 );
  bool wasforced=false;
  for(unsigned i=0; i<getNumberOfVessels(); ++i) {
    if( (functions[i]->applyForce( tmpforces )) ) {
      wasforced=true;
      for(unsigned j=0; j<forcesToApply.size(); ++j) forcesToApply[j]+=tmpforces[j];
    }
  }
  return wasforced;
}

void ActionWithVessel::retrieveDomain( std::string& min, std::string& max ) {
  plumed_merror("If your function is periodic you need to add a retrieveDomain function so that ActionWithVessel can retrieve the domain");
}

Vessel* ActionWithVessel::getVesselWithName( const std::string& mynam ) {
  int target=-1;
  for(unsigned i=0; i<functions.size(); ++i) {
    if( functions[i]->getName().find(mynam)!=std::string::npos ) {
      if( target<0 ) target=i;
      else error("found more than one " + mynam + " object in action");
    }
  }
  return functions[target].get();
}

}
}
