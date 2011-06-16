#ifdef __PLUMED_HAS_MATHEVAL
#include "ActionRegister.h"
#include "Function.h"
#include <cassert>
#include <matheval.h>

using namespace std;

namespace PLMD{


//+PLUMEDOC FUNCTION MATHEVAL
/**
Combination of more CVs using a matheval expression
   
* VAR tells the names of the variables used in the FUNC string
  If absent, the names are set to x,y,z. With more than 3 arguments explicit
  names are compulsory
* FUNC gives the function definition in matheval syntax

\verbatim
MATHEVAL LABEL=c1 ARG=distance1,distance2 VAR=a,b FUNC=(a+b)/10.0
MATHEVAL LABEL=c2 ARG=distance1,distance2         FUNC=(x-y)/sin(y)
\endverbatim


*/
//+ENDPLUMEDOC


class FunctionMatheval :
  public Function
{
  void* evaluator;
  vector<void*> evaluator_deriv;
  vector<string> var;
  string func;
  vector<double> values;
  vector<char*> names;
public:
  FunctionMatheval(const ActionOptions&);
  ~FunctionMatheval();
  void calculate();
};


PLUMED_REGISTER_ACTION(FunctionMatheval,"MATHEVAL")

FunctionMatheval::FunctionMatheval(const ActionOptions&ao):
Action(ao),
Function(ao),
evaluator_deriv(getNumberOfArguments()),
values(getNumberOfArguments()),
names(getNumberOfArguments())
{
  parseVector("VAR",var);
  if(var.size()==0){
    var.resize(getNumberOfArguments());
    assert(getNumberOfArguments()<=3);
    if(var.size()>0) var[0]="x";
    if(var.size()>1) var[1]="y";
    if(var.size()>2) var[2]="z";
  }
  assert(var.size()==getNumberOfArguments());
  parse("FUNC",func);
  checkRead();

  evaluator=evaluator_create(const_cast<char*>(func.c_str()));

  char **check_names;
  int    check_count;
  evaluator_get_variables(evaluator,&check_names,&check_count);
  assert(check_count==int(getNumberOfArguments()));
  for(unsigned i=0;i<getNumberOfArguments();i++){
    bool found=false;
    for(unsigned j=0;j<getNumberOfArguments();j++){
      if(var[i]==check_names[j])found=true;
    }
    assert(found);
  }

  for(unsigned i=0;i<getNumberOfArguments();i++)
    evaluator_deriv[i]=evaluator_derivative(evaluator,const_cast<char*>(var[i].c_str()));

  addValueWithDerivatives("");

  log.printf("  with function : %s\n",func.c_str());
  log.printf("  with variables :");
  for(unsigned i=0;i<var.size();i++) log.printf(" %s",var[i].c_str());
  log.printf("\n");
}

void FunctionMatheval::calculate(){
  for(unsigned i=0;i<getNumberOfArguments();i++) values[i]=getArgument(i);
  for(unsigned i=0;i<getNumberOfArguments();i++) names[i]=const_cast<char*>(var[i].c_str());
  setValue(evaluator_evaluate(evaluator,names.size(),&names[0],&values[0]));

  for(unsigned i=0;i<getNumberOfArguments();i++){
    setDerivatives(i,evaluator_evaluate(evaluator_deriv[i],names.size(),&names[0],&values[0]));
  }
}

FunctionMatheval::~FunctionMatheval(){
  evaluator_destroy(evaluator);
  for(unsigned i=0;i<evaluator_deriv.size();i++)evaluator_destroy(evaluator_deriv[i]);
}

}

#endif

