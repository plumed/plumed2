#include "ColvarModifier.h"
#include "ColvarModifierFunctions.h"

namespace PLMD {

//+PLUMEDOC MODIFIER MORE_THAN
/**

Calculate the number of colvars that are greater than the value specified.  To make this quantity continuous it is calculated using:

\f[
\sigma = 1.0 - \sum_i \frac{ 1 - \left(\frac{ s_i }{ r_0 }\right)^{nn} }{ 1 - \left(\frac{ s_i }{ r_0 }\right)^{mm} }
\f]

where \f$r_0\f$ is the value specied in the MORE_THAN keyword.  By default \f$nn\f$ and \f$mm\f$ are set equal to 6 and 12 respectively.
These values can be adjusted by using the LOGIC_NN and LOGIC_MM keywords.  Once calculated the final value is stored inside a value called
<label>.less_than<\f$r_0\f$>.

*/
//+ENDPLUMEDOC


ColvarModifierMore::ColvarModifierMore(ColvarWithModifiers* act):
ColvarModifier(act),
total(0)
{
  double r_0; parse("MORE_THAN", false, r_0);
  int nn=6; parse("LOGIC_NN", true, nn);
  int mm=12; parse("LOGIC_MM", true, mm);
  sf.set( nn, mm, r_0, 0.0 );
  std::string target; Tools::convert( r_0, target );
  std::string snn; Tools::convert( nn, snn );
  std::string smm; Tools::convert( mm, smm );
  report("  number of values more than " + target + ".  Switching function parameters are " + snn + " " + smm ); 
  addValue("less_than" + target); 
}

void ColvarModifierMore::finishCalculation(){
  setValue(0, total, 1.0 );
  total=0;
}

double ColvarModifierMore::differential( const unsigned& ival, const double& value ){
  assert(ival==0); double df;
  total+=( 1.0 - sf.calculate(value, df) );
  return -df*value;
} 

}
