/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "ActionWithValue.h"
#include "ActionWithArguments.h"
#include "tools/Stopwatch.h"
#include "tools/Exception.h"
#include "tools/OpenMP.h"

using namespace std;
namespace PLMD {

void ActionWithValue::registerKeywords(Keywords& keys) {
  keys.setComponentsIntroduction("By default the value of the calculated quantity can be referenced elsewhere in the "
                                 "input file by using the label of the action.  Alternatively this Action can be used "
                                 "to calculate the following quantities by employing the keywords listed "
                                 "below.  These quanties can be referenced elsewhere in the input by using this Action's "
                                 "label followed by a dot and the name of the quantity required from the list below.");
  keys.addFlag("NUMERICAL_DERIVATIVES", false, "calculate the derivatives for these quantities numerically");
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize");
  keys.addFlag("TIMINGS",false,"output information on the timings of the various parts of the calculation");
}

void ActionWithValue::noAnalyticalDerivatives(Keywords& keys) {
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addFlag("NUMERICAL_DERIVATIVES",true,"analytical derivatives are not implemented for this keyword so numerical derivatives are always used");
}

void ActionWithValue::componentsAreNotOptional(Keywords& keys) {
  keys.setComponentsIntroduction("By default this Action calculates the following quantities. These quanties can "
                                 "be referenced elsewhere in the input by using this Action's label followed by a "
                                 "dot and the name of the quantity required from the list below.");
}

void ActionWithValue::useCustomisableComponents(Keywords& keys) {
  keys.setComponentsIntroduction("The names of the components in this action can be customized by the user in the "
                                 "actions input file.  However, in addition to these customizable components the "
                                 "following quantities will always be output");
}

ActionWithValue::ActionWithValue(const ActionOptions&ao):
  Action(ao),
  noderiv(true),
  numericalDerivatives(false),
  serial(false),
  timers(false),
  nactive_tasks(0)
{
  if( keywords.exists("NUMERICAL_DERIVATIVES") ) parseFlag("NUMERICAL_DERIVATIVES",numericalDerivatives);
  if(numericalDerivatives) log.printf("  using numerical derivatives\n");
  parseFlag("SERIAL",serial); parseFlag("TIMINGS",timers);
  stopwatch.start(); stopwatch.pause();
}

ActionWithValue::~ActionWithValue() {
  stopwatch.start(); stopwatch.stop();
  if(timers) {
    log.printf("timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log<<stopwatch;
  }
}

const std::string & ActionWithValue::getLabelOfActionThatCalculates() const {
  const ActionWithArguments* aa = dynamic_cast<const ActionWithArguments*>( this );
  if( aa ){
      if( !aa->done_over_stream ) return getLabel(); 
  }
  return getLabel();
}

void ActionWithValue::addActionToChain( ActionWithValue* act ){
  for(const auto & p : actions_to_do_after ){ 
      if( act->getLabel()==p->getLabel() ) return;
  }
  actions_to_do_after.push_back( act );
}

void ActionWithValue::clearInputForces() {
  for(unsigned i=0; i<values.size(); i++) values[i]->clearInputForce();
}

void ActionWithValue::clearDerivatives( const bool& force ) {
  const ActionWithArguments* aa = dynamic_cast<const ActionWithArguments*>( this );
  if( aa ){ if( !force && aa->done_over_stream ) return; }
  unsigned nt = OpenMP::getGoodNumThreads(values);
  #pragma omp parallel num_threads(nt)
  {
    #pragma omp for
    for(unsigned i=0; i<values.size(); i++) values[i]->clearDerivatives();
  }
  for(const auto & p : actions_to_do_after ) p->clearDerivatives( true ); 
}

// -- These are the routine for copying the value pointers to other classes -- //

bool ActionWithValue::exists( const std::string& name ) const {
  for(unsigned i=0; i<values.size(); ++i) {
    if (values[i]->name==name) return true;
  }
  return false;
}

Value* ActionWithValue::copyOutput( const std::string& name ) const {
  for(unsigned i=0; i<values.size(); ++i) {
    if (values[i]->name==name) return values[i].get();
  }
  plumed_merror("there is no pointer with name " + name);
  return NULL;
}

Value* ActionWithValue::copyOutput( const unsigned& n ) const {
  plumed_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n].get();
}

// -- HERE WE HAVE THE STUFF FOR THE DEFAULT VALUE -- //

void ActionWithValue::addValue( const std::vector<unsigned>& shape ) {
  plumed_massert(values.empty(),"You have already added the default value for this action");
  values.emplace_back(new Value(this,getLabel(), false, shape ) );
}

void ActionWithValue::addValueWithDerivatives( const std::vector<unsigned>& shape ) {
  if( shape.size()>0 && shape.size()!=getNumberOfDerivatives() ) plumed_merror("should not be adding non zero rank values with derivatives");
  plumed_massert(values.empty(),"You have already added the default value for this action");
  values.emplace_back(new Value(this,getLabel(), true, shape ) );
}

void ActionWithValue::setNotPeriodic() {
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->min=0; values[0]->max=0;
  values[0]->setupPeriodicity();
}

void ActionWithValue::setPeriodic( const std::string& min, const std::string& max ) {
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->setDomain( min, max );
}

Value* ActionWithValue::getPntrToValue() {
  plumed_dbg_massert(values.size()==1,"The number of components is not equal to one");
  plumed_dbg_massert(values[0]->name==getLabel(), "The value you are trying to retrieve is not the default");
  return values[0].get();
}

// -- HERE WE HAVE THE STUFF FOR NAMED VALUES / COMPONENTS -- //

void ActionWithValue::addComponent( const std::string& name, const std::vector<unsigned>& shape ) {
  if( !keywords.outputComponentExists(name,true) ) {
    warning("a description of component " + name + " has not been added to the manual. Components should be registered like keywords in "
            "registerKeywords as described in the developer docs.");
  }
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0; i<values.size(); ++i) {
    plumed_massert(values[i]->name!=getLabel(),"Cannot mix single values with components");
    plumed_massert(values[i]->name!=thename&&name!="bias","Since PLUMED 2.3 the component 'bias' is automatically added to all biases by the general constructor!\n"
                   "Remove the line addComponent(\"bias\") from your bias.");
    plumed_massert(values[i]->name!=thename,"there is already a value with this name");
  }
  values.emplace_back(new Value(this,thename, false, shape ) );
  std::string msg="  added component to this action:  "+thename+" \n";
  log.printf(msg.c_str());
}

void ActionWithValue::addComponentWithDerivatives( const std::string& name, const std::vector<unsigned>& shape ) {
  if( !keywords.outputComponentExists(name,true) ) {
    warning("a description of component " + name + " has not been added to the manual. Components should be registered like keywords in "
            "registerKeywords as described in the developer doc.");
  }
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0; i<values.size(); ++i) {
    plumed_massert(values[i]->name!=getLabel(),"Cannot mix single values with components");
    plumed_massert(values[i]->name!=thename&&name!="bias","Since PLUMED 2.3 the component 'bias' is automatically added to all biases by the general constructor!\n"
                   "Remove the line addComponentWithDerivatives(\"bias\") from your bias.");
    plumed_massert(values[i]->name!=thename,"there is already a value with this name");
  }
  values.emplace_back(new Value(this,thename, true, shape ) );
  std::string msg="  added component to this action:  "+thename+" \n";
  log.printf(msg.c_str());
}

int ActionWithValue::getComponent( const std::string& name ) const {
  plumed_massert( !exists( getLabel() ), "You should not be calling this routine if you are using a value");
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0; i<values.size(); ++i) {
    if (values[i]->name==thename) return i;
  }
  plumed_merror("there is no component with name " + name);
  return -1;
}

std::string ActionWithValue::getComponentsList( ) const {
  std::string complist;
  for(unsigned i=0; i<values.size(); ++i) {
    complist+=values[i]->name+" ";
  }
  return complist;
}

std::vector<std::string> ActionWithValue::getComponentsVector( ) const {
  std::vector<std::string> complist;
  for(unsigned i=0; i<values.size(); ++i) {
    complist.push_back(values[i]->name);
  }
  return complist;
}

void ActionWithValue::componentIsNotPeriodic( const std::string& name ) {
  int kk=getComponent(name);
  values[kk]->min=0; values[kk]->max=0;
  values[kk]->setupPeriodicity();
}

void ActionWithValue::componentIsPeriodic( const std::string& name, const std::string& min, const std::string& max ) {
  int kk=getComponent(name);
  values[kk]->setDomain(min,max);
}

void ActionWithValue::setGradientsIfNeeded() {
  if(isOptionOn("GRADIENTS")) {
    for(unsigned i=0; i<values.size(); i++) values[i]->setGradients();
  }
}

void ActionWithValue::turnOnDerivatives() {
  // Turn on the derivatives
  noderiv=false;
  // Resize the derivatives
  for(unsigned i=0; i<values.size(); ++i) values[i]->resizeDerivatives( getNumberOfDerivatives() );
  // And turn on the derivatives in all actions on which we are dependent
  for(unsigned i=0; i<getDependencies().size(); ++i) {
    ActionWithValue* vv=dynamic_cast<ActionWithValue*>( getDependencies()[i] );
    if(vv) vv->turnOnDerivatives();
  }
}

Value* ActionWithValue::getPntrToOutput( const unsigned& ind ) const {
  plumed_dbg_massert(ind<values.size(),"you have requested a pointer that is out of bounds");
  return values[ind];
}

Value* ActionWithValue::getPntrToComponent( const std::string& name ){
  int kk=getComponent(name);
  return values[kk].get();
}

Value* ActionWithValue::getPntrToComponent( int n ) {
  plumed_dbg_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n].get();
}

void ActionWithValue::interpretDataLabel( const std::string& mystr, ActionWithArguments* myuser, std::vector<Value*>& args ){
  // Check for streams
  unsigned nstr=0; for(unsigned i=0;i<values.size();++i){ if( values[i]->getRank()>0 ) nstr++; }

  if( mystr=="" || mystr==getLabel() ){
      // Retrieve value with name of action
      if( values[0]->name!=getLabel() ) myuser->error("action " + getLabel() + " does not have a value");
      args.push_back( values[0] ); values[0]->interpretDataRequest( myuser->getLabel(), "" ); 
  } else if( mystr==getLabel() + ".*" ){
      // Retrieve all scalar values
      retrieveAllScalarValuesInLoop( args );
      for(unsigned j=0;j<args.size();++j) args[j]->interpretDataRequest( myuser->getLabel(), "" );
  } else if( mystr.find(".")!=std::string::npos && exists( mystr ) ) {
      // Retrieve value with specific name
      args.push_back( copyOutput( mystr ) ); copyOutput( mystr )->interpretDataRequest( myuser->getLabel(), "" ); 
  } else {
      std::size_t dot1 = mystr.find_first_of('.'); std::string thelab = mystr.substr(0,dot1); 
      plumed_assert( thelab==getLabel() ); std::string rest = mystr.substr(dot1+1); 
      if( rest.find_first_of('.')==std::string::npos ){
         plumed_assert( values.size()==1 ); plumed_assert( values[0]->getRank()>0 && values[0]->getName()==getLabel() );
         args.push_back( values[0] ); values[0]->interpretDataRequest( myuser->getLabel(), rest ); 
      } else {
         std::size_t dot2 = rest.find_first_of('.'); std::string thecomp = rest.substr(0,dot2);
         if( !exists( thelab + "." + thecomp ) ) myuser->error("could not find component with label " + thelab + "." + thecomp );
         args.push_back( copyOutput( thelab + "." + thecomp ) ); 
         getPntrToComponent( thecomp )->interpretDataRequest( myuser->getLabel(), rest.substr(dot2+1) );
      }
  }
}

void ActionWithValue::addTaskToList( const unsigned& taskCode ) {
  fullTaskList.push_back( taskCode ); taskFlags.push_back(0);
  plumed_assert( fullTaskList.size()==taskFlags.size() );
}

void ActionWithValue::selectActiveTasks( std::vector<unsigned>& tflags ){
  buildCurrentTaskList( tflags );
  for(const auto & p : actions_to_do_after ) p->selectActiveTasks( tflags );
//  for(unsigned i=0;i<values.size();++i){
//      if( values[i]->getRank()>0 ) values[i]->activateTasks( tflags );
//  }
}

void ActionWithValue::runAllTasks() {
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) { stride=1; rank=0; }

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nactive_tasks ) nt=nactive_tasks/stride/10;
  if( nt==0 || !threadSafe() ) nt=1;

  // Build the list of active tasks 
  taskFlags.assign(taskFlags.size(),0); selectActiveTasks( taskFlags ); nactive_tasks = 0;
  for(unsigned i=0; i<fullTaskList.size(); ++i) {
    if( taskFlags[i]>0 ) nactive_tasks++;
  }
  unsigned n=0; partialTaskList.resize( nactive_tasks ); indexOfTaskInFullList.resize( nactive_tasks );
  for(unsigned i=0; i<fullTaskList.size(); ++i) {
    // Deactivate sets inactive tasks to number not equal to zero
    if( taskFlags[i]>0 ) {
      partialTaskList[n] = fullTaskList[i]; indexOfTaskInFullList[n]=i; n++;
    }
  }

  // Get the total number of streamed quantities that we need
  unsigned nquantities = 0; getNumberOfStreamedQuantities( nquantities );
  for(const auto & p : actions_to_do_after ) p->getNumberOfStreamedQuantities( nquantities );
  // Get size for buffer
  unsigned bufsize=0; getSizeOfBuffer( nactive_tasks, bufsize );
  if( buffer.size()!=bufsize ) buffer.resize( bufsize );
  // Clear buffer
  buffer.assign( buffer.size(), 0.0 );

  // Recover the number of derivatives we require
  unsigned nderivatives = 0;
  if( !noderiv ){
      nderivatives = getNumberOfDerivatives();
      for(const auto & p : actions_to_do_after ){
          unsigned nnd = p->getNumberOfDerivatives();
          if( nnd>nderivatives ) nderivatives = nnd;
      }
  }

  if(timers) stopwatch.start("2 Loop over tasks");
  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_buffer;
    if( nt>1 ) omp_buffer.resize( bufsize, 0.0 );
    MultiValue myvals( nquantities, nderivatives  );
    myvals.clearAll();

    #pragma omp for nowait
    for(unsigned i=rank; i<nactive_tasks; i+=stride) {
      // Calculate the stuff in the loop for this action
      runTask( indexOfTaskInFullList[i], partialTaskList[i], myvals );

      // Now transfer the data to the actions that accumulate values from the calculated quantities
      if( nt>1 ) {
        gatherAccumulators( indexOfTaskInFullList[i], myvals, omp_buffer );
      } else {
        gatherAccumulators( indexOfTaskInFullList[i], myvals, buffer );
      }

      // Clear the value
      myvals.clearAll();
    }
    #pragma omp critical
    if(nt>1) for(unsigned i=0; i<bufsize; ++i) buffer[i]+=omp_buffer[i];
  }
  if(timers) stopwatch.stop("2 Loop over tasks");

  if(timers) stopwatch.start("3 MPI gather");
  // MPI Gather everything
  if( !serial && buffer.size()>0 ) comm.Sum( buffer );
  // Update the elements that are makign contributions to the sum here
  // this causes problems if we do it in prepare
  if(timers) stopwatch.stop("3 MPI gather");

  if(timers) stopwatch.start("4 Finishing computations");
  finishComputations( buffer );
  if(timers) stopwatch.stop("4 Finishing computations");
}

unsigned ActionWithValue::getNumberOfStreamedQuantities( unsigned& nquants ) const {
  for(unsigned i=0;i<values.size();++i){ values[i]->streampos=nquants; nquants++; }
  for(const auto & p : actions_to_do_after ) p->getNumberOfStreamedQuantities( nquants );
  return nquants;
}

void ActionWithValue::getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize ){
  for(unsigned i=0;i<values.size();++i){
      values[i]->bufstart=bufsize; 
      if( values[i]->getRank()==0 && !noderiv ) bufsize += 1 + values[i]->getNumberOfDerivatives();
      else if( values[i]->getRank()==0 ) bufsize += 1;
      else if( values[i]->storedata ){
          if( values[i]->hasDeriv ) bufsize += values[i]->getSize(); else bufsize += nactive_tasks;
      }
      // if( values[i]->getRank()==0 ) bufsize += 1 + values[i]->getNumberOfDerivatives();
  }
  for(const auto & p : actions_to_do_after ) p->getSizeOfBuffer( nactive_tasks, bufsize );
}

void ActionWithValue::runTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  myvals.setTaskIndex(task_index); performTask( current, myvals ); 
  for(const auto & p : actions_to_do_after ){
     if( p->isActive() ) p->runTask( task_index, current, myvals );
  }
}

void ActionWithValue::rerunTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  for(const auto & p : getDependencies() ){
     ActionWithArguments* aa=dynamic_cast<ActionWithArguments*>(p);
     if( aa ){
        ActionWithValue* av=dynamic_cast<ActionWithValue*>(p);
        av->rerunTask( task_index, current, myvals );
     }
  }
  myvals.setTaskIndex(task_index); performTask( current, myvals );
}

void ActionWithValue::gatherAccumulators( const unsigned& taskCode, const MultiValue& myvals, std::vector<double>& buffer ) const {
  for(unsigned i=0;i<values.size();++i){
      unsigned sind = values[i]->streampos, bufstart = values[i]->bufstart; 
      if( values[i]->getRank()==0 ){
           buffer[bufstart] += myvals.get(sind);
           if( values[i]->hasDerivatives() ){
               for(unsigned k=0;k<myvals.getNumberActive();++k){
                   unsigned kindex = myvals.getActiveIndex(k); buffer[bufstart + 1 + kindex] += myvals.getDerivative(sind,kindex);
               }
           }
      } else if( values[i]->storedata ){
           unsigned vindex = bufstart + taskCode*(1+values[i]->getNumberOfDerivatives()); buffer[vindex] += myvals.get(sind);
      } 
  }
  for(const auto & p : actions_to_do_after ){
     if( p->isActive() ) p->gatherAccumulators( taskCode, myvals, buffer );
  }
}

void ActionWithValue::retrieveAllScalarValuesInLoop( std::vector<Value*>& myvals ){
  for(unsigned i=0;i<values.size();++i){
      if( values[i]->getRank()==0 ){
          bool found=false;
          for(unsigned j=0;j<myvals.size();++j){
              if( values[i]->getName()==myvals[j]->getName() ){ found=true; break; }
          }
          if( !found ) myvals.push_back( values[i] );
      }
  }
  for(const auto & p : actions_to_do_after ) p->retrieveAllScalarValuesInLoop( myvals );
}

void ActionWithValue::finishComputations( const std::vector<double>& buffer ){
  for(unsigned i=0;i<values.size();++i){
      unsigned bufstart = values[i]->bufstart; 
      if( values[i]->reset ) values[i]->data.assign( values[i]->data.size(), 0 );
      if( values[i]->storedata ){
          for(unsigned j=0;j<values[i]->getSize();++j) values[i]->add( j, buffer[bufstart+j] ); 
      }
      if( !doNotCalculateDerivatives() && values[i]->getRank()==0 ){ 
          for(unsigned j=0;j<values[i]->getNumberOfDerivatives();++j) values[i]->setDerivative( j, buffer[bufstart+1+j] );
      }
  } 
  for(const auto & p : actions_to_do_after ){
     if( p->isActive() ) p->finishComputations( buffer );
  }
}

}
