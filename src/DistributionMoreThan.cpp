#include "DistributionFunctions.h"

namespace PLMD {

void more_than::writeDocs(std::string& docs){
  std::ostringstream ostr;
  ostr<<"\\par MORE_THAN\n\n";
  ostr<<"Calculate the number of functions that are more than a target value.  To make this quantity continuous it is calculated using:\n";
  ostr<<"\\f$\n";
  ostr<<"\\sigma = 1.0 - \\sum_i \\frac{ 1 - \\left(\\frac{ s_i }{ r_0 }\\right)^{n} }{ 1 - \\left(\\frac{ s_i }{ r_0 }\\right)^{m} }\n";
  ostr<<"\\f$\n";
  ostr<<", where \\f$r_0\\f$ is specied after the LESS_THAN keyword.  By default \\f$n\\f$ and \\f$m\\f$ are set equal to 6 and 12 respectively.\n";
  ostr<<"These values can be adjusted by using using three arguments for less than (\\f$r_0\\f$, \\f$n\\f$, \\f$m\\f$).  Once calculated the final value is referenced\n";
  ostr<<"using label.lt\\f$r_0\\f$.\n"; 
  docs=ostr.str();
}

more_than::more_than( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  if( parameters.size()==3 ){
     Tools::convert(parameters[0],r_0);
     Tools::convert(parameters[1],nn); Tools::convert(parameters[2],mm);
     sf.set(nn, mm, r_0, 0.0 );
  } else if(parameters.size()==1){
     Tools::convert(parameters[0],r_0);
     nn=6; mm=12;
     sf.set(nn, mm, r_0, 0.0 );
  } else {
     error("MORE_THAN keyword takes one of three arguments");
  }
}

std::string more_than::message(){
  std::ostringstream ostr;
  ostr<<"number of values more than "<<r_0<<".  Switching function parameters are "<<nn<<" and "<<mm;
  return ostr.str();
}

void more_than::calculate( Value* value_in, std::vector<Value>& aux, Value* value_out ){
  copyDerivatives( value_in, value_out );
  double p, df, f; p=value_in->get(); 
  f=1.0 - sf.calculate(p, df); df*=-p;
  value_out->chainRule(df); value_out->set(f);
}

void more_than::finish( const double& p, Value* value_out ){
  value_out->set(p);
}

}
