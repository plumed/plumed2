#include "DistributionFunctions.h"

namespace PLMD {

std::string less_than::documentation(){
  std::ostringstream ostr;
  ostr<<SwitchingFunction::documentation();
  return ostr.str();
}

less_than::less_than( const std::string& parameters ) :
DistributionFunction(parameters)
{
  std::string errormsg;
  sf.set( parameters, errormsg ); 
  if( errormsg.size()!=0 ) error( errormsg ); 
  addAccumulator( true );
}

std::string less_than::message(){
  std::ostringstream ostr;
  ostr<<"number of values less than "<<sf.description();
  return ostr.str();
}

void less_than::printKeywords( Log& log ){
  sf.printKeywords( log );
}

std::string less_than::getLabel(){
  std::string vv;
  Tools::convert( sf.get_r0(), vv );
  return "lt" + vv;
}

void less_than::calculate( Value* value_in, std::vector<Value>& aux ){
  copyValue( 0, value_in );
  double p, df, f; p=value_in->get(); 
  f=sf.calculate(p, df); df*=p;
  chainRule( 0, df ); setValue( 0, f );
}

void less_than::finish( Value* value_out ){
  extractDerivatives( 0, value_out );
  value_out->set( getPntrToAccumulator(0)->get() );
}

}
