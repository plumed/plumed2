#include "ActionWithDistribution.h"

using namespace std;
using namespace PLMD;

void ActionWithDistribution::registerKeywords(Keywords& keys){
  keys.add("optional","NL_STRIDE","the frequency with which the neighbor list should be updated.");
  keys.add("optional","NL_TOL","when accumulating sums quantities that contribute less than this will be ignored.");
}

void ActionWithDistribution::autoParallelize(Keywords& keys){
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize over collective variables");
}

ActionWithDistribution::ActionWithDistribution(const ActionOptions&ao):
  Action(ao),
  read(false),
  all_values(true),
  serial(true),
  updateFreq(0),
  lastUpdate(0),
  reduceAtNextStep(false),
  tolerance(0)
{
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  if(serial)log.printf("  doing calculation in serial\n");
  if( keywords.exists("NL_STRIDE") ) parse("NL_STRIDE",updateFreq);
  if( keywords.exists("NL_TOL") ){
      if(updateFreq>0 && keywords.exists("NL_TOL") ){ tolerance=epsilon; parse("NL_TOL",tolerance); }
      log.printf("  updating calculation every %d steps.  Ignoring contributions less than %lf\n",updateFreq,tolerance);
  }
}

ActionWithDistribution::~ActionWithDistribution(){
  for(unsigned i=0;i<functions.size();++i) delete functions[i];
}

void ActionWithDistribution::addDistributionFunction( std::string name, DistributionFunction* fun ){
  if(all_values) all_values=false;  // Possibly will add functionality to delete all values here

  // Check function is good 
  if( !fun->check() ) error( fun->errorMessage() );

  // Add a value
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  a->addComponentWithDerivatives( name );
  unsigned fno=a->getNumberOfComponents()-1;
  final_values.push_back( a->copyOutput( fno ) );

  // Add the function   
  plumed_massert( fno==functions.size(), "Number of functions does not match number of values" );
  functions.push_back( fun );
  log.printf("  value %s contains %s\n",( (a->copyOutput( fno ))->getName() ).c_str(),( functions[fno]->message() ).c_str() );
}

void ActionWithDistribution::requestDistribution(){
  read=true; bool dothis; std::vector<std::string> params;
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  plumed_massert(a,"can only do distribution on ActionsWithValue");

  // We prepare the first step as if we are doing a neighbor list update
  prepareForNeighborListUpdate(); reduceAtNextStep=true;

  if(all_values){
     error("No function has been specified");
  } else {
     plumed_massert( functions.size()==final_values.size(), "number of functions does not match number of values" );
     // This sets up the dynamic list that holds what we are calculating
     for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){ members.addIndexToList(i); }
     members.activateAll(); members.updateActiveMembers();
  }
}

void ActionWithDistribution::prepare(){
 if(reduceAtNextStep){
    completeNeighborListUpdate();
    // Setup the functions by declaring enough space to hold the derivatives
    for(unsigned i=0;i<functions.size();++i) functions[i]->setNumberOfDerivatives( final_values[0]->getNumberOfDerivatives() );
    // Setup the buffers for mpi gather
    if(!serial){
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
    for(unsigned i=0;i<functions.size();++i) functions[i]->setNumberOfDerivatives( final_values[0]->getNumberOfDerivatives() );
    // Setup the buffers for mpi gather
    if(!serial){
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
  // Create a value to store stuff in 
  Value* tmpvalue=new Value();

  unsigned kk; double contribution; bool keep;
  for(unsigned i=rank;i<members.getNumberActive();i+=stride){
      // Retrieve the function we are calculating from the dynamic list
      kk=members[i];
      // Make sure we have enough derivatives in this value
      unsigned nder=getThisFunctionsNumberOfDerivatives(kk);
      if( tmpvalue->getNumberOfDerivatives()!=nder ) tmpvalue->resizeDerivatives( nder );
 
      // Calculate the value of this particular function 
      calculateThisFunction( kk, tmpvalue, aux );

      // Skip if we are not calculating this particular value
      if( reduceAtNextStep && !tmpvalue->valueHasBeenSet() ){ members.deactivate(kk); deactivate(kk); continue; }

      // Now incorporate the derivative of the function into the derivatives for the min etc
      if( updateFreq>0 ){ keep=false; } else { keep=true; }
      for(unsigned j=0;j<functions.size();++j){
         functions[j]->calculate( tmpvalue, aux );
         if( functions[j]->sizableContribution( tolerance ) ){ 
             keep=true; functions[j]->mergeDerivatives( kk, *this );
         }
      }
      tmpvalue->clearDerivatives();
      // If the contribution of this quantity is very small at neighbour list time ignore it
      // untill next neighbour list time
      if( reduceAtNextStep && !keep ){ members.deactivate(kk); deactivate(kk); } 
  }
  // Update the dynamic list 
  if(reduceAtNextStep){ members.mpi_gatherActiveMembers( comm ); }
  // MPI Gather everything
  if(!serial){ 
     unsigned bufsize=0;
     for(unsigned i=0;i<functions.size();++i) functions[i]->copyDataToBuffers( bufsize, buffer );
     plumed_assert( bufsize==buffer.size() ); 
     comm.Sum( &buffer[0],buffer.size() ); 
     bufsize=0;
     for(unsigned i=0;i<functions.size();++i) functions[i]->retrieveDataFromBuffers( bufsize, buffer );
     plumed_assert( bufsize==buffer.size() );
  }

  // Delete the tmpvalues
  delete tmpvalue; 
  // Set the final value of the function
  for(unsigned j=0;j<final_values.size();++j) functions[j]->finish( final_values[j] ); 
}
