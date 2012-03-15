#include "DistributionFunctions.h"

namespace PLMD {

void less_than::writeDocs(std::string& docs){
  std::ostringstream ostr;
  ostr<<"\\par LESS_THAN\n\n"; 
  ostr<<"Calculate the number of functions that are less than a target value using a switching function. "<<SwitchingFunction::documentation()<<"\n";
  docs=ostr.str();
}

less_than::less_than( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  printf("Hello from less than %d\n",parameters.size() );
  if( parameters.size()==1 ){
     sf.set( parameters[0] ); 
  } else {
     error("LESS_THAN keyword takes one switching function as argument");
  }
  addAccumulator( true );
}

std::string less_than::message(){
  std::ostringstream ostr;
  ostr<<"number of values less than "<<sf.description();
  return ostr.str();
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
