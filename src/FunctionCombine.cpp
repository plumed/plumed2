#include "ActionRegister.h"
#include "Function.h"

#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC FUNCTION COMBINE
/**
Calculate a polynomial combination of a set of other variables

\par Examples
The following input tells plumed to print the distance between atoms 3 and 5
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
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(FunctionCombine,"COMBINE")

void FunctionCombine::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.add("compulsory","COEFFICIENTS","1.0","the coefficients of the arguments in your function");
  keys.add("compulsory","POWERS","1.0","the powers to which you are raising each of the arguments in your function");
  keys.addFlag("NORMALIZE",false,"normalize all the coefficents so that in total they are equal to one");
}

FunctionCombine::FunctionCombine(const ActionOptions&ao):
Action(ao),
Function(ao),
normalize(false),
coefficients(getNumberOfArguments(),1.0),
powers(getNumberOfArguments(),1.0)
{
  printf("In function combine\n");
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
 
  addValueWithDerivatives("");
  double min(0),max(0); std::vector<std::string> period;
  parseVector("PERIODIC",period);
  if(period.size()==1 && period[0]=="NO"){
    getValue("")->setPeriodicity(false);
  } else if(period.size()==2 && Tools::convert(period[0],min) && Tools::convert(period[1],max)){
    getValue("")->setPeriodicity(true);
    getValue("")->setDomain(min,max);
  } else error("missing PERIODIC keyword");
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
  setValue(combine);
}

}


