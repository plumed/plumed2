#include "DistributionFunctions.h"

namespace PLMD {

//+PLUMEDOC MODIFIER MIN
/**

Calculate the minimum from the set of defined colvars.  To make this quantity continuous the minimum is calculated using:

\f[
\textrm{min} = \frac{\beta}{ \log \sum_i \exp\left( \frac{\beta}{s_i} \right) }
\f]

The keyword MIN=\f$\beta\f$ tells plumed to calculate this quantity and sets \f$\beta\f$ to an appropriate value.  Once 
calcualted the final value is referneced using label.min  

*/
//+ENDPLUMEDOC

min::min( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  if( parameters.size()==1 ){
     Tools::convert(parameters[0],beta);
  } else {
     error("MIN keyword should have one argument - the value for beta");
  }
}

std::string min::message(){
  std::ostringstream ostr;
  ostr<<"the minimum value. beta is equal to "<<beta; 
  return ostr.str();
}

double min::calculate( Value* value_in, std::vector<Value>& aux, Value* value_out ){
  copyDerivatives( value_in, value_out );
  double p, df, tmp; p=value_in->get();
  tmp=exp( beta/p ); df=tmp/(p*p); 
  value_out->chainRule(df); value_out->set(tmp);
  return tmp;
}

void min::finish( const double& p, Value* value_out ){
  double df, dist=beta/std::log(p); df=dist*dist/p;
  value_out->chainRule(df); value_out->set(dist);
}

}
