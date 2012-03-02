#include "DistributionFunctions.h"

namespace PLMD {

void sum::writeDocs(std::string& docs){
  std::ostringstream ostr;
  ostr<<"\\par SUM \n\n";
  ostr<<"Calculate the sum of all the colvars in the distribution.  Once calculated the final value is referenced\n";
  ostr<<"using lable.sum.\n";
  docs=ostr.str();
}

sum::sum( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  plumed_massert(parameters.size()==1,"should pass one parameter total values");
  Tools::convert(parameters[0],prev_nval);
}

std::string sum::message(){
  std::ostringstream ostr;
  ostr<<"the sum of all the values"; 
  return ostr.str();
}

void sum::calculate( Value* value_in, std::vector<Value>& aux, Value* value_out ){
  copyDerivatives( value_in, value_out ); value_out->set( value_in->get() );
  nval+=1.0;
}

void sum::gather( PlumedCommunicator& comm ){
  comm.Sum(&nval,1); 
}

void sum::finish( const double& p, Value* value_out ){
  if ( nval!=prev_nval ) printf("WARNING: A neighbor list is causing discontinuities in a sum");
  prev_nval=nval; nval=0;
  value_out->set(p);
}

}
