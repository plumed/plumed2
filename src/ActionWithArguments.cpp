#include "ActionWithArguments.h"
#include "ActionWithValue.h"
#include "PlumedMain.h"
#include "ActionSet.h"

using namespace std;
using namespace PLMD;

void ActionWithArguments::parseArgumentList(const std::string&key,std::vector<Value*>&arg){
  vector<string> c;
  arg.clear();
  parseVector(key,c);
  for(unsigned i=0;i<c.size();i++){
    std::size_t dot=c[i].find_first_of('.');
    if(dot!=string::npos){
      string a=c[i].substr(0,dot);
      string name=c[i].substr(dot+1);
      if(a=="*"){
        plumed_massert(name=="*","arguments in the form *.something are not allowed, but for *.*");
        std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
        for(unsigned j=0;j<all.size();j++){
          for(int k=0;k<all[j]->getNumberOfValues();++k){
            arg.push_back(all[j]->getValue(k));
          }
        };
      } else {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(a);
        plumed_massert(action,"cannot find action named " + a);
        if(name=="*"){
          vector<string> s=action->getValueNames();
          for(unsigned j=0;j<s.size();j++)arg.push_back(action->getValue(s[j]));
        } else {
          plumed_massert(action->hasNamedValue(name),"action " + a + " has no component with name "+ name);
          arg.push_back(action->getValue(name));
        }
      }
    } else if(c[i]=="*"){
      std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
      for(unsigned j=0;j<all.size();j++){
        plumed_massert(all[j]->hasNamedValue(""),"action " + all[j]->getLabel() + " has no default component (unnamed one)");
        arg.push_back(all[j]->getValue(""));
      };
    } else{
      ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(c[i]);
      plumed_massert(action,"cannot find action named " +c[i]);
      plumed_massert(action->hasNamedValue(""),"action "+c[i]+" has no default component (unnamed one)");
      arg.push_back(action->getValue(""));
    }
  }
}

void ActionWithArguments::requestArguments(const vector<Value*> &arg){
  plumed_massert(!lockRequestArguments,"requested argument list can only be changed in the prepare() method");
  arguments=arg;
  clearDependencies();
  for(unsigned i=0;i<arguments.size();i++) addDependency(&arguments[i]->getAction());
}

ActionWithArguments::ActionWithArguments(const ActionOptions&ao):
  Action(ao),
  lockRequestArguments(false)
{
  vector<Value*> arg;
  parseArgumentList("ARG",arg);

  if(arg.size()>0){
    log.printf("  with arguments");
    for(unsigned i=0;i<arg.size();i++) log.printf(" %s",arg[i]->getFullName().c_str());
    log.printf("\n");
  }

  requestArguments(arg);

}

void ActionWithArguments::calculateNumericalDerivatives(){
  ActionWithValue*a=dynamic_cast<ActionWithValue*>(this);
  plumed_massert(a,"cannot compute numerical derivatives for an action without values");
  const int nval=a->getNumberOfValues();
  const int npar=a->getNumberOfParameters();
  std::vector<double> value (nval*npar);
  for(int i=0;i<npar;i++){
    double arg0=arguments[i]->get();
    arguments[i]->set(arg0+sqrt(epsilon));
    calculate();
    arguments[i]->set(arg0);
    for(int j=0;j<nval;j++){
      value[i*nval+j]=a->getValue(j)->get();
    }
  }
  calculate();
  std::vector<double> value0(nval);
  for(int j=0;j<nval;j++){
    Value* v=a->getValue(j);
    if(v->hasDerivatives())for(int i=0;i<npar;i++) v->setDerivatives(i,(value[i*nval+j]-a->getValue(j)->get())/sqrt(epsilon));
  }
}






