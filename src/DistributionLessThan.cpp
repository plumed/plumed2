#include "DistributionFunctions.h"

namespace PLMD {

//+PLUMEDOC MODIFIER LESS_THAN
/**

Calculate the number of functions that are less than a target value.  To make this quantity continuous it is calculated using:
\f$
\sigma = \sum_i \frac{ 1 - \left(\frac{ s_i }{ r_0 }\right)^{n} }{ 1 - \left(\frac{ s_i }{ r_0 }\right)^{m} }
\f$
, where \f$r_0\f$ is specied after the LESS_THAN keyword.  By default \f$n\f$ and \f$m\f$ are set equal to 6 and 12 respectively.
These values can be adjusted by using using three arguments for less than (\f$r_0\f$, \f$n\f$, \f$m\f$).  Once calculated the final value is referenced
using label.lt\f$r_0\f$. 

*/
//+ENDPLUMEDOC


less_than::less_than( const std::vector<std::string>& parameters ) :
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
     error("LESS_THAN keyword takes one of three arguments");
  }
}

std::string less_than::message(){
  std::ostringstream ostr;
  ostr<<"number of values less than "<<r_0<<".  Switching function parameters are "<<nn<<" and "<<mm;
  return ostr.str();
}

double less_than::compute( const double p, double& df ){
  double f; f=sf.calculate(p, df); df*=p;
  return f;
}

double less_than::last_step( const double p, double& df ){
  df=1.0; return p;
}

}
