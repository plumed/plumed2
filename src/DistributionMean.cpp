#include "DistributionFunctions.h"

namespace PLMD {

void mean::writeDocs(std::string& docs){
  std::ostringstream ostr;
  ostr<<"\\par AVERAGE \n\n";
  ostr<<"Calculate the average of all the colvars in the distribution.  Once calculated the final value is referenced\n";
  ostr<<"using lable.average.\n";
  docs=ostr.str();
}

mean::mean( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  plumed_massert(parameters.size()==1,"should pass one parameter total values");
  Tools::convert(parameters[0],prev_nval); nval=0.;
}

std::string mean::message(){
  std::ostringstream ostr;
  ostr<<"the average value"; 
  return ostr.str();
}

void mean::calculate( Value* value_in, std::vector<Value>& aux, Value* value_out ){
  copyDerivatives( value_in, value_out ); value_out->set( value_in->get() );
  nval+=1.0;
}

void mean::gather( PlumedCommunicator& comm ){
  comm.Sum(&nval,1); 
}

void mean::finish( const double& p, Value* value_out ){
  if ( nval!=prev_nval ) printf("WARNING: A neighbor list is causing discontinuities in an average");
  prev_nval=nval; 
  value_out->chainRule(1.0/nval); value_out->set(p/nval); 
  nval=0;
}

}
