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
//  keys.add("optional","HISTOGRAM", "create a discretized histogram of the distribution");
}

ActionWithDistribution::ActionWithDistribution(const ActionOptions&ao):
  Action(ao),
  read(false),
  all_values(true)
{
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
     totals.resize( a->getNumberOfComponents() );
     plumed_massert( functions.size()==final_values.size(), "number of functions does not match number of values" );
  }
}

void ActionWithDistribution::calculate(){
  plumed_massert( read, "you must have a call to readDistributionKeywords somewhere" );  

  if(all_values){
      for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i) calculateThisFunction( i, final_values[i] );
  } else {
      // Reset all totals
      for(unsigned j=0;j<totals.size();++j) totals[j]=0.0;
      // Create a value to store stuff in 
      Value* tmpvalue=new Value();
      Value* tmp2value=new Value();
      unsigned kk;
      for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){
          // Make sure we have enough derivatives in this value
          unsigned nder=getThisFunctionsNumberOfDerivatives(i);
          if( tmpvalue->getNumberOfDerivatives()!=nder ){
              tmpvalue->resizeDerivatives( nder );
              tmp2value->resizeDerivatives( nder );
          }
 
          // Calculate the value of this particular function 
          calculateThisFunction( i, tmpvalue );
          // Now incorporate the derivative of the function into the derivatives for the min etc
          for(unsigned j=0;j<totals.size();++j){
             totals[j]+=functions[j]->calculate( tmpvalue, tmp2value );
             mergeDerivatives( i, tmp2value, final_values[j] );
             tmp2value->clearDerivatives();
          }
          tmpvalue->clearDerivatives();
      }
      //
      delete tmpvalue; delete tmp2value;
      // Set the final value of the function
      for(unsigned j=0;j<totals.size();++j) functions[j]->finish( totals[j], final_values[j] );
  }
}
