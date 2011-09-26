#include "ActionRegister.h"
#include "Function.h"

#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC FUNCTION COMBINE
/**
Calculate the polynomial combination of other variables

\par Syntax
\verbatim
COMBINE ARG=x1,x2,... [POWERS=p1,p2,...] [COEFFICIENTS=c1,c2,...] [NORMALIZE]
\endverbatim
The resulting variable has value
\f$
  \sum_i c_i x_i^{p_i}
\f$. When not present, powers and coefficient are implicitly equal to 1.
If NORMALIZE is present, the c coefficients are first normalized.


\par Example
The following input is printing the distance between atoms 3 and 5
its square (as computed from the x,y,z components) and the distance
again as computed from the square root of the square.
\verbatim
DISTANCE LABEL=dist      ATOMS=3,5 COMPONENTS
COMBINE  LABEL=distance2 ARG=dist.x,dist.y,dist.z POWERS=2,2,2
COMBINE  LABEL=distance  ARG=distance2 POWERS=0.5
PRINT ARG=distance,distance2
\endverbatim
(See also \ref PRINT and \ref DISTANCE).


*/
//+ENDPLUMEDOC


class FunctionCombine :
  public Function
{
  bool normalize;
  std::vector<double> coefficients;
  std::vector<double> powers;
public:
  FunctionCombine(const ActionOptions&);
  void calculate();
};


PLUMED_REGISTER_ACTION(FunctionCombine,"COMBINE")

FunctionCombine::FunctionCombine(const ActionOptions&ao):
Function(ao),
normalize(false),
coefficients(getNumberOfArguments(),1.0),
powers(getNumberOfArguments(),1.0)
{
  registerKeyword(1,"COEFFICIENTS","the coefficients for the terms in the sum");
  registerKeyword(1,"POWERS","the powers to raise each colvar value to");
  registerKeyword(0,"NORMALIZE","normalize the coefficients");
  readFunction();

  parseVector("COEFFICIENTS",coefficients);
  assert(coefficients.size()==static_cast<unsigned>(getNumberOfArguments()));
  parseVector("POWERS",powers);
  assert(powers.size()==static_cast<unsigned>(getNumberOfArguments()));

  parseFlag("NORMALIZE",normalize);

  if(normalize){
    double n=0.0;
    for(unsigned i=0;i<coefficients.size();i++) n+=coefficients[i];
    for(unsigned i=0;i<coefficients.size();i++) coefficients[i]*=(1.0/n);
  }

  //addValueWithDerivatives("");

  checkRead();

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
  //setValue(combine);
}

}


