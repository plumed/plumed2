#include "DistributionFunctions.h"

namespace PLMD {

mean::mean( const std::string& parameters ) :
DistributionFunction(parameters)
{
  Tools::convert(parameters,nval); 
  addAccumulator( true );
  addAccumulator( false );
}

std::string mean::message(){
  std::ostringstream ostr;
  ostr<<"the average value"; 
  return ostr.str();
}

void mean::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in ); 
  setValue( 1, 1.0 );
}

void mean::finish( Value* value_out ){
  if ( getPntrToAccumulator(1)->get()!=nval ) printf("WARNING: A neighbor list is causing discontinuities in an average");
  extractDerivatives( 0, value_out );
  value_out->chainRule(1.0/nval); value_out->set(getPntrToAccumulator(0)->get()/nval); 
}

}
