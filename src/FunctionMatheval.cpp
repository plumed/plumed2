#include "ActionRegister.h"
#include "Function.h"
#include <cassert>

#ifdef __PLUMED_HAS_MATHEVAL
#include <matheval.h>
#endif

using namespace std;

namespace PLMD{


//+PLUMEDOC FUNCTION MATHEVAL
/*
Calculate a combination of variables using a matheval expression.

\par Examples
The following input tells plumed to print the angle between vectors
identified by atoms 1,2 and atoms 2,3
its square (as computed from the x,y,z components) and the distance
again as computed from the square root of the square.
\verbatim
DISTANCE LABEL=d1 ATOMS=1,2 COMPONENTS
DISTANCE LABEL=d2 ATOMS=2,3 COMPONENTS
MATHEVAL ...
  LABEL=theta
  ARG=d1.x,d1.y,d1.z,d2.x,d2.y,d2.z
  VAR=ax,ay,az,bx,by,bz
  FUNC=acos((ax*bx+ay*by+az*bz)/sqrt((ax*ax+ay*ay+az*az)*(bx*bx+by*by+bz*bz))
... MATHEVAL
PRINT ARG=theta
\endverbatim
(See also \ref PRINT and \ref DISTANCE).

\attention
The MATHEVAL object only works if libmatheval is installed on the system and
PLUMED has been linked to it

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
  static void registerKeywords(Keywords& keys);
};

#ifdef __PLUMED_HAS_MATHEVAL
PLUMED_REGISTER_ACTION(FunctionMatheval,"MATHEVAL")

void FunctionMatheval::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG"); keys.use("PERIODIC");
  keys.add("compulsory","FUNC","the function you wish to evaluate");
  keys.add("optional","VAR","the names to give each of the arguments in the function.  If you have up to three arguments in your function you can use x, y and z to refer to them.  Otherwise you must use this flag to give your variables names.");
}

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
  addValueWithDerivatives(); 
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
    setDerivative(i,evaluator_evaluate(evaluator_deriv[i],names.size(),&names[0],&values[0]));
  }
}

FunctionMatheval::~FunctionMatheval(){
  evaluator_destroy(evaluator);
  for(unsigned i=0;i<evaluator_deriv.size();i++)evaluator_destroy(evaluator_deriv[i]);
}

#endif

}


