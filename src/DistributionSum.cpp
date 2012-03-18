#include "DistributionFunctions.h"

namespace PLMD {

sum::sum( const std::string& parameters ) :
DistributionFunction(parameters)
{
  Tools::convert(parameters,nval); 
  addAccumulator( true );
  addAccumulator( false );
}

std::string sum::message(){
  std::ostringstream ostr;
  ostr<<"the sum of all the values"; 
  return ostr.str();
}

void sum::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in ); 
  setValue( 1, 1.0 );
}

void sum::finish( Value* value_out ){
  if ( getPntrToAccumulator(1)->get()!=nval ) printf("WARNING: A neighbor list is causing discontinuities in a sum");
  extractDerivatives( 0, value_out ); value_out->set( getPntrToAccumulator(0)->get() );
}

}
