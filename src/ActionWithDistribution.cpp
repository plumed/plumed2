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
     if( a->checkNumericalDerivatives() && getNumberOfFunctionsInDistribution()>1 ){
         error("Can only use numerical derivatives for distribution functions or single colvars");
     }
     std::string ss; 
     for(int i=0;i<getNumberOfFunctionsInDistribution();++i){
        Tools::convert(i,ss);
        a->addComponentWithDerivatives( "value" + ss );
        final_values.push_back( a->copyOutput(i) );
     }
  } else {
     plumed_massert( functions.size()==final_values.size(), "number of functions does not match number of values" );
     totals.resize( a->getNumberOfComponents() );
     // This sets up the dynamic list that holds what we are calculating
     for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){ members.addIndexToList(i); }
     members.activateAll(); members.updateActiveMembers();
  }
}

void ActionWithDistribution::prepare(){
 if(reduceAtNextStep){
    completeNeighborListUpdate();
    reduceAtNextStep=false;
 }
 if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
    members.activateAll();
    members.updateActiveMembers();
    prepareForNeighborListUpdate();
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
  if(all_values){
      // It is not at all straightforward to parallelize this.  Also is it worth it?
      for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){ calculateThisFunction( i, final_values[i], aux ); }
  } else {
      // Reset all totals
      for(unsigned j=0;j<totals.size();++j) totals[j]=0.0;
      // Create a value to store stuff in 
      Value* tmpvalue=new Value();
      Value* tmp2value=new Value();

      unsigned kk; double contribution; bool keep;
      for(unsigned i=rank;i<members.getNumberActive();i+=stride){
          // Retrieve the function we are calculating from the dynamic list
          kk=members[i];
          // Make sure we have enough derivatives in this value
          unsigned nder=getThisFunctionsNumberOfDerivatives(kk);
          if( tmpvalue->getNumberOfDerivatives()!=nder ){
              tmpvalue->resizeDerivatives( nder );
              tmp2value->resizeDerivatives( nder );
          }
 
          // Calculate the value of this particular function 
          calculateThisFunction( kk, tmpvalue, aux );

          // Skip if we are not calculating this particular value
          if( reduceAtNextStep && !tmpvalue->valueHasBeenSet() ){ members.deactivate(kk); deactivate(kk); continue; }

          // Now incorporate the derivative of the function into the derivatives for the min etc
          keep=false;
          for(unsigned j=0;j<totals.size();++j){
             contribution=functions[j]->calculate( tmpvalue, aux, tmp2value );
             if( updateFreq>0 || contribution>=tolerance ){ keep=true; totals[j]+=contribution; } 
             mergeDerivatives( kk, tmp2value, final_values[j] );
             tmp2value->clearDerivatives();
          }
          tmpvalue->clearDerivatives();
          // If the contribution of this quantity is very small at neighbour list time ignore it
          // untill next neighbour list time
          if( reduceAtNextStep && !keep ){ members.deactivate(kk); deactivate(kk); } 
      }
      // Update the dynamic list 
      if(reduceAtNextStep){ members.mpi_gatherActiveMembers( comm ); }
      // MPI Gather the totals and derivatives
      if(!serial){ 
         comm.Sum( &totals[0],totals.size() ); 
         for(unsigned j=0;j<totals.size();++j){ final_values[j]->gatherDerivatives( comm ); }
      }

      // Delete the tmpvalues
      delete tmpvalue; delete tmp2value;
      // Set the final value of the function
      for(unsigned j=0;j<totals.size();++j) functions[j]->finish( totals[j], final_values[j] ); 
  }

  for(unsigned i=0;i<totals.size();++i){
      if( !final_values[i]->valueHasBeenSet() ) error("error in input for this action.  Not all values have been set"); 
  }

}
