#include "ColvarModifier.h"
#include "ColvarModifierFunctions.h"

namespace PLMD {

//+PLUMEDOC MODIFIER MEAN 
/**

Calculate the average for the distribution.  Once calculated the average is stored inside a value called <label>.average.

*/
//+ENDPLUMEDOC


ColvarModifierMean::ColvarModifierMean(ColvarWithModifiers* act):
ColvarModifier(act),
total(0)
{
  parse("NCOLVAR",false,ncv); assert(ncv>0);
  report("  using average value"); 
  addValue("average");
}

void ColvarModifierMean::finishCalculation(){
  setValue(0, total/ncv, 1.0/ncv);
  total=0;
}

double ColvarModifierMean::differential( const unsigned& ival, const double& value ){
  assert(ival==0); 
  total+=value;
  return 1.0;
} 

}
