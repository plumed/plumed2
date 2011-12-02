#include "ColvarModifier.h"
#include "ColvarModifierFunctions.h"

namespace PLMD {

//+PLUMEDOC MODIFIER SUM
/**

Calculate the sum of all the colvars in the distribution.  Once calculated this sum is stored inside a value called <lable>.sum.

*/
//+ENDPLUMEDOC


ColvarModifierSum::ColvarModifierSum(ColvarWithModifiers* act):
ColvarModifier(act),
total(0)
{
  report("  using sum of all values");
  addValue("sum");
}

void ColvarModifierSum::finishCalculation(){
  setValue(0, total, 1.0);
  total=0;
}

double ColvarModifierSum::differential( const unsigned& ival, const double& value ){
  assert(ival==0); 
  total+=value;
  return 1.0;
} 

}
