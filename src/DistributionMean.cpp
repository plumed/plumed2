#include "DistributionFunctions.h"

namespace PLMD {

//+PLUMEDOC MODIFIER AVERAGE
/**

Calculate the average of all the colvars in the distribution.  Once calculated the final value is referenced
using lable.average.

*/
//+ENDPLUMEDOC


mean::mean( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  plumed_massert(parameters.size()==1,"should pass one parameter total values");
  Tools::convert(parameters[0],nvalues);
}

std::string mean::message(){
  std::ostringstream ostr;
  ostr<<"the average value"; 
  return ostr.str();
}

double mean::calculate( Value* value_in, std::vector<Value>& aux, Value* value_out ){
  copyDerivatives( value_in, value_out ); value_out->set( value_in->get() );
  return value_in->get();
}

void mean::finish( const double& p, Value* value_out ){
  value_out->chainRule(1.0/nvalues); value_out->set(p/nvalues); 
}

}
