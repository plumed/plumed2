#include "ActionWithArguments.h"
#include "ActionWithValue.h"
#include "PlumedMain.h"

using namespace std;
using namespace PLMD;

ActionWithArguments::ActionWithArguments(const ActionOptions&ao):
  ActionWithValue(ao),
  lockRequestArguments(false)
{
  registerKeyword(2, "ARG", "a list of plumed actions that provide the input to this action");
  allowKeyword("ARG");
}

void ActionWithArguments::readActionWithArguments( const std::vector<double>& domain ){
  vector<string> c; arguments.clear();
  parseVector("ARG",c);
  for(unsigned i=0;i<c.size();i++){
    std::size_t dot=c[i].find_first_of('.');
    if(dot!=string::npos){
      string a=c[i].substr(0,dot);
      string name=c[i].substr(dot+1);
      if(a=="*"){
        assert(name=="*");
        std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
        for(unsigned j=0;j<all.size();j++){
          for(int k=0;k<all[j]->getNumberOfValues();++k){
            arguments.push_back( all[j]->getValuePointer(k) );
          }
        };
      } else {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(a);
        assert(action);
        if(name=="*"){
          // vector<string> s=action->getValueNames();
          for(unsigned j=0;j<action->getNumberOfValues();j++) arguments.push_back( action->getValuePointer(j) );
        } else {
          // assert(action->hasNamedValue(name));
          arguments.push_back( action->getValuePointer(name) );
        }
      }
    } else if(c[i]=="*"){
      std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
      for(unsigned j=0;j<all.size();j++){
        assert(all[j]->getNumberOfValues()==1); 
        arguments.push_back( all[j]->getValuePointer(0) );
      }
    } else{
      ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(c[i]);
      assert(action); assert(action->getNumberOfValues()==1); 
      arguments.push_back( action->getValuePointer(0) );
    }
  }

  if(arguments.size()>0){
    clearDependencies();
    log.printf("  with arguments");
    for(unsigned i=0;i<arguments.size();i++){
        log.printf(" %s",arguments[i]->myname.c_str());
        addDependency(&arguments[i]->getAction());
    }
    log.printf("\n");
  }
  readActionWithValue( arguments.size(), domain );
}

void ActionWithArguments::calculateNumericalDerivatives(){
  const int nval=getNumberOfValues();
  const int npar=arguments.size(); 
  std::vector<double> value (nval*npar);
  for(int i=0;i<npar;i++){
    double arg0=getArgument(i); 
    arguments[i]->set( arg0+sqrt(epsilon) );
    calculate();
    arguments[i]->set( arg0 );
    for(int j=0;j<nval;j++){
      value[i*nval+j]=getValue(j);
    }
  }
  calculate();
  clearDerivatives();
  //std::vector<double> value0(nval);
  for(int j=0;j<nval;j++){
    //Value* v=a->getValue(j);
    if ( isMeaningfulToDifferentiate(j) ){
       for(int i=0;i<npar;i++) addDerivative( j, i, ( value[i*nval+j] - getValue(j) )/sqrt(epsilon) );
    }
  }
}

void ActionWithArguments::printArgumentNames( FILE* fp ){
  for(unsigned i=0;i<arguments.size();i++){
    fprintf(fp," %s",arguments[i]->myname.c_str());
  }
}






