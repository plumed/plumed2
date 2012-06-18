#include "ActionRegister.h"
#include "Function.h"

#include <cmath>
#include <cassert>
#include <algorithm>
#include <utility>

using namespace std;

namespace PLMD{

//+PLUMEDOC FUNCTION SORT
/**
Just sorts its arguments

\par Examples
The following input tells plumed to print the distance of the closest and of
the farthest atoms to atom 1, chosen among atoms from 2 to 5
\verbatim
d12:  DISTANCE ATOMS=1,2
d13:  DISTANCE ATOMS=1,3
d14:  DISTANCE ATOMS=1,4
d15:  DISTANCE ATOMS=1,5
sort: SORT ARG=d12,d13,d14,d15
PRINT ARG=sort.1,sort.4
\endverbatim
(See also \ref PRINT and \ref DISTANCE).

*/
//+ENDPLUMEDOC


class FunctionSort :
  public Function
{
public:
  FunctionSort(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(FunctionSort,"SORT")

void FunctionSort::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
}

FunctionSort::FunctionSort(const ActionOptions&ao):
Action(ao),
Function(ao)
{
  for(unsigned i=0;i<getNumberOfArguments();++i){
    string s;
    Tools::convert(i+1,s);
    if(getPntrToArgument(i)->isPeriodic())
      error("Cannot sort periodic values (check argument "+s+")");
    addComponentWithDerivatives(s); 
    getPntrToComponent(i)->setNotPeriodic();
  }
  checkRead();

}

void FunctionSort::calculate(){
  vector<pair<double,int> > vals(getNumberOfArguments());
  for(unsigned i=0;i<getNumberOfArguments();++i){
    vals[i].first=getArgument(i);
// In this manner I remember from which argument the component depends:
    vals[i].second=i;
  }
// STL sort sorts based on first element (value) then second (index)
  sort(vals.begin(),vals.end());
  for(unsigned i=0;i<getNumberOfComponents();++i){
    Value* v=getPntrToComponent(i);
    v->set(vals[i].first);
    setDerivative(v,vals[i].second,1.0);
  }
}

}


