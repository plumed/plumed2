#include "Function.h"
#include "Colvar.h"

using namespace PLMD;
using namespace std;

Function::Function(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
ActionWithArguments(ao)
{
    setNumberOfParameters(getNumberOfArguments());
    double min(0),max(0); std::vector<std::string> period;
    parseVector("PERIODIC",period);
    if(period.size()==1 && period[0]=="NO"){
      getValue("")->setPeriodicity(false);
    } else if(period.size()==2 && Tools::convert(period[0],min) && Tools::convert(period[1],max)){
      getValue("")->setPeriodicity(true);
      getValue("")->setDomain(min,max);
    } else error("missing PERIODIC keyword");
    checkRead();
}

void Function::apply(){

  vector<double>   f(getNumberOfArguments(),0.0);
  bool at_least_one_forced=false;

  for(int i=0;i<getNumberOfValues();++i){
    if(!getValue(i)->checkForced()) continue;
    at_least_one_forced=true;
    const vector<double> & derivatives(getValue(i)->getDerivatives());
    for(unsigned j=0;j<derivatives.size();j++){
      f[j]+=getForce(i)*derivatives[j];
    }
  }
  if(at_least_one_forced) for(unsigned i=0;i<getNumberOfArguments();++i){
    getArguments()[i]->addForce(f[i]);
  }
}

void Function::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.add("compulsory","PERIODIC","nosize","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
}


