#include "DistributionFunctions.h"

namespace PLMD {

std::string within::documentation(){
  std::ostringstream ostr;
  ostr<<"To make this quantity continuous it is calculated using "<<HistogramBead::documentation(false);
  return ostr.str();
}

within::within( const std::string& parameters ) :
DistributionFunction(parameters)
{ 
  hist.set( parameters, "" );
  addAccumulator( true );
}

std::string within::message(){
  std::ostringstream ostr;
  ostr<<"number of values "<<hist.description(); 
  return ostr.str();
}

void within::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in ); 
  double df, f; f=hist.calculate( value_in->get() , df );
  chainRule(0, df); setValue(0, f);
}

void within::finish( Value* value_out ){
  extractDerivatives( 0, value_out );
  value_out->set( getPntrToAccumulator(0)->get() );
}

}
