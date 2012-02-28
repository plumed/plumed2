#include "ActionWithDistribution.h"

using namespace std;
using namespace PLMD;

void ActionWithDistribution::registerKeywords(Keywords& keys){
  keys.add("nohtml","MIN","take the minimum value from these variables.  The continuous version of the minimum described above is calculated and beta must be specified in input");
//  keys.add("optional","MAX", "take the maximum value from these variables");
  keys.add("nohtml","SUM", "take the sum of these variables");
  keys.add("nohtml","AVERAGE", "take the average value of these variables");
  keys.add("nohtml","LESS_THAN", "take the number of variables less than the specified target.  This quantity is made differentiable using a switching function.  You can control the parameters of this switching function by specifying three numbers to the keyword (r_0, nn and mm).  If you are happy with the default values of nn=6 and mm=12 then you need only specify the target r_0.  The number of values less than the target is stored in a value called lt<target>.");
  keys.add("nohtml","MORE_THAN", "take the number of variables more than the specified target.  This quantity is made differentiable using a switching function.  You can control the parameters of this switching function by specifying three numbers to the keyword (r_0, nn and mm).  If you are happy with the default values of nn=6 and mm=12 then you need only specify the target r_0.  The number of values less than the target is stored in a value called gt<target>.");
  keys.add("nohtml","HISTOGRAM", "create a discretized histogram of the distribution.  This keyword's input should consist of one or two numbers");
  keys.add("nohtml","RANGE", "the range in which to calculate the histogram");
  keys.add("nohtml", "WITHIN", "The number of values within a certain range.  This keyword's input should consist of one or two numbers.");
  keys.reserve("optional","NL_STRIDE","the frequency with which the neighbour list should be updated");
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
  reduceAtNextStep(false)
{
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  if(serial)log.printf("  doing calculation in serial\n");
  if( keywords.exists("NL_STRIDE") ) parse("NL_STRIDE",updateFreq);
  if(updateFreq>0) reduceAtNextStep=true;
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

void ActionWithDistribution::readDistributionKeywords(){
  read=true; bool dothis; std::vector<std::string> params;
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  plumed_massert(a,"can only do distribution on ActionsWithValue");

  // Read SUM keyword
  parseFlag("SUM",dothis);
  if( dothis ){
     addDistributionFunction( "sum", new sum(params) );   
  }
  // Read AVERAGE keyword
  parseFlag("AVERAGE",dothis);
  if( dothis ){
     params.resize(1); Tools::convert(getNumberOfFunctionsInDistribution(),params[0]);
     addDistributionFunction( "average", new mean(params) );
  }
  // Read MIN keyword
  parseVector("MIN",params);
  if( params.size()!=0 ){ 
     addDistributionFunction( "min", new min(params) );
  }
  // Read Less_THAN keyword
  parseVector("LESS_THAN",params);
  if( params.size()!=0 ){
     addDistributionFunction( "lt" + params[0], new less_than(params) );
  }
  // Read MORE_THAN keyword
  parseVector("MORE_THAN",params);
  if( params.size()!=0 ){
     addDistributionFunction( "gt" + params[0], new more_than(params) );
  }
  // Read HISTOGRAM keyword
  parseVector("HISTOGRAM",params);
  if( params.size()!=0 ){
      std::vector<double> range(2); parseVector("RANGE",range);
      if(range[1]<range[0]) error("range is meaningless");
      int nbins; Tools::convert(params[0],nbins);
      std::vector<std::string> hparams(2);
      if(params.size()==2){
          hparams.resize(3);
          hparams[2]=params[1];
      } else if(params.size()!=1){
          error("Histogram keyword should either specify just the number"
                " of bins or the number of bins and the ammount of smearing"); 
      }
      double lb,ub,delr=(range[1]-range[0])/static_cast<double>(nbins);
      for(int i=0;i<nbins;++i){
          lb=range[0]+i*delr; Tools::convert( lb, hparams[0] );
          ub=range[0]+(i+1)*delr; Tools::convert( ub, hparams[1] );
          addDistributionFunction( "between" + hparams[0] + "&" +hparams[1], new within(hparams) );
      }
  }
  // Read within keywords
  parseVector("WITHIN",params);
  if( params.size()!=0 ){
      addDistributionFunction( "between" + params[0] + "&" +params[1], new within(params) );
  } else {
      for(unsigned i=1;;++i){
         if( !parseNumberedVector("WITHIN",i,params) ) break;
         addDistributionFunction( "between" + params[0] + "&" +params[1], new within(params) );
      }
  }

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

//void ActionWithDistribution::resetMembers(){
//  members.activateAll(); 
//  members.updateActiveMembers();
//}

void ActionWithDistribution::prepare(){
 if(reduceAtNextStep){
    completeNeighbourListUpdate();
    reduceAtNextStep=false;
 }
 if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
    members.activateAll();
    members.updateActiveMembers();
    prepareForNeighbourListUpdate();
    reduceAtNextStep=true;
    lastUpdate=getStep();
 }
}

void ActionWithDistribution::calculate(){
  plumed_massert( read, "you must have a call to readDistributionKeywords somewhere" );  

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

      unsigned kk;
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
          if( reduceAtNextStep && !tmpvalue->valueHasBeenSet() ){ members.deactivate(kk); continue; }
          // Now incorporate the derivative of the function into the derivatives for the min etc
          for(unsigned j=0;j<totals.size();++j){
             totals[j]+=functions[j]->calculate( tmpvalue, aux, tmp2value );
             mergeDerivatives( kk, tmp2value, final_values[j] );
             tmp2value->clearDerivatives();
          }
          tmpvalue->clearDerivatives();
      }
      // Update the dynamic list during neighbour list update
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
