#include "ActionRegister.h"
#include "Function.h"

#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC FUNCTION COMBINE
/**
Combinatin of more CVs
   
\par Syntax


*/
//+ENDPLUMEDOC


class FunctionCombine :
  public Function
{
  std::vector<double> coefficients;
  std::vector<double> powers;
public:
  FunctionCombine(const ActionOptions&);
  void calculate();
};


PLUMED_REGISTER_ACTION(FunctionCombine,"COMBINE")

FunctionCombine::FunctionCombine(const ActionOptions&ao):
Action(ao),
Function(ao),
coefficients(getNumberOfArguments(),1.0),
powers(getNumberOfArguments(),1.0)
{
  parseVector("COEFFICIENTS",coefficients);
  assert(coefficients.size()==static_cast<unsigned>(getNumberOfArguments()));
  parseVector("POWERS",powers);
  assert(powers.size()==static_cast<unsigned>(getNumberOfArguments()));
  checkRead();

  addValueWithDerivatives("");
  log.printf("  with coefficients:");
  for(unsigned i=0;i<coefficients.size();i++) log.printf(" %f",coefficients[i]);
  log.printf("\n");
  log.printf("  and powers:");
  for(unsigned i=0;i<powers.size();i++) log.printf(" %f",powers[i]);
  log.printf("\n");
}

void FunctionCombine::calculate(){
  double combine=0.0;
  for(unsigned i=0;i<coefficients.size();++i){
    combine+=coefficients[i]*pow(getArgument(i),powers[i]);
    setDerivatives(i,coefficients[i]*powers[i]*pow(getArgument(i),powers[i]-1.0));
  };
  setValue(combine);
}

}


