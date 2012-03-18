#include "DistributionFunctions.h"

namespace PLMD {

std::string more_than::documentation(){
  std::ostringstream ostr;
  ostr<<"This is calculated using \\f$1.0 - s(r)\\f$, where \\f$s(r)\\f$ is a switching function. ";
  ostr<<SwitchingFunction::documentation();
  return ostr.str();
}

more_than::more_than( const std::string& parameters ) :
DistributionFunction(parameters)
{
  sf.set( parameters );
  addAccumulator( true );
}

std::string more_than::message(){
  std::ostringstream ostr;
  ostr<<"number of values more than "<<sf.description();
  return ostr.str();
}

void more_than::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in );
  double p, df, f; p=value_in->get(); 
  f=1.0 - sf.calculate(p, df); df*=-p;
  chainRule( 0, df ); setValue( 0, f );
}

void more_than::finish( Value* value_out ){
  extractDerivatives( 0, value_out );
  value_out->set( getPntrToAccumulator(0)->get() );
}

}
