#include "ActionRegister.h"
#include "Function.h"

#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC FUNCTION COMBINE
/**
Calculate a polynomial combination of other variables - i.e. \f$ \sum_i c_i x_i^{p_i}

\f$. When not present, powers and coefficient are implicitly equal to 1.
If NORMALIZE is present, the c coefficients are first normalized.

\par Example
The following computes two distances along with a linear comination of the two.
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=5,7 LABEL=d2
COMBINE LABLE=c1 ARG=d1,d2 POWERS=1,1 COEFFICIENTS=0.5,0.75
\endverbatim

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
  registerKeyword(0,"COEFFICIENTS","(default=1) the coefficients for the terms in the sum");
  registerKeyword(0,"POWERS","(default=1) the powers to raise each colvar value to");
  registerKeyword(0,"NORMALIZE","normalize the coefficients");
  readFunction();

  parseVector("COEFFICIENTS",coefficients);
  if(coefficients.size()==0){
      coefficients.resize( getNumberOfArguments() );
      for(unsigned i=0;i<coefficients.size();++i) coefficients[i]=1.0;
  } else if (coefficients.size()!=static_cast<unsigned>(getNumberOfArguments())){
      error("number of coefficents does not match number of arguments");
  }
  parseVector("POWERS",powers);
  if(powers.size()==0){
      powers.resize( getNumberOfArguments() );
      for(unsigned i=0;i<powers.size();++i) powers[i]=1.0;
  } else if (powers.size()!=static_cast<unsigned>(getNumberOfArguments())){
      error("number of powers does not match number of arguments");
  }
  parseFlag("NORMALIZE",normalize);

  if(normalize){
    double n=0.0;
    for(unsigned i=0;i<coefficients.size();i++) n+=coefficients[i];
    for(unsigned i=0;i<coefficients.size();i++) coefficients[i]*=(1.0/n);
  }
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
}

}


