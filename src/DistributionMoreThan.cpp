#include "DistributionFunctions.h"

namespace PLMD {

void more_than::writeDocs(std::string& docs){
  std::ostringstream ostr;
  ostr<<"\\par MORE_THAN\n\n";
  ostr<<"Calculate the number of functions that are more than a target value using \\f$1.0 - s(r)\\f$, where \\f$s(r)\\f$ is a switching function. ";
  ostr<<SwitchingFunction::documentation()<<"\n";
  docs=ostr.str();
}

more_than::more_than( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  if( parameters.size()==1 ){
     sf.set( parameters[0] );
  } else {
     error("MORE_THAN keyword takes one switching function as argument");
  }
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
