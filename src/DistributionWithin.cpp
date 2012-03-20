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
  std::string errormsg;
  hist.set( parameters, "", errormsg );
  if( errormsg.size()!=0 ) error( errormsg );
  addAccumulator( true );
}

std::string within::message(){
  std::ostringstream ostr;
  ostr<<"number of values "<<hist.description(); 
  return ostr.str();
}

void within::printKeywords( Log& log ){
  hist.printKeywords( log );
}

std::string within::getLabel(){
  std::string lb, ub;
  Tools::convert( hist.getlowb(), lb );
  Tools::convert( hist.getbigb(), ub );
  return "between" + lb + "&" + ub;
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
