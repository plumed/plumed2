#include "ColvarModifierFunctions.h"

namespace PLMD {

//+PLUMEDOC MODIFIER MIN
/**

Calculate the minimum from the set of defined colvars.  To make this quantity continuous the minimum is calculated using:

\f[
\textrm{min} = \frac{\beta}{ \log \sum_i \exp\left( \frac{\beta}{s_i} \right) }
\f]

By default \f$\beta\f$ is set equal to 50 but you can change this quantity by using the BETA keyword.  Once calculated the 
minimum is stored inside a value called <label>.min.

*/
//+ENDPLUMEDOC


ColvarModifierMin::ColvarModifierMin(ColvarWithModifiers* act):
ColvarModifier(act),
total(0),
beta(50)
{
  parse("BETA",true,beta);
  std::string sbeta; Tools::convert(beta,sbeta);
  report("  using minimum value.  Value of beta is " + sbeta );
  addValue("min");
}

void ColvarModifierMin::finishCalculation(){
  double dist=beta/std::log(total); 
  setValue(0, dist, dist*dist/total );
  total=0;
}

double ColvarModifierMin::differential( const unsigned& ival, const double& value ){
  assert(ival==0);
  double tmp=exp( beta/value ); total+=tmp; 
  return tmp/(value*value);
} 

}
