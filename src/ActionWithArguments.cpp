#include "ActionWithArguments.h"
#include "PlumedMain.h"

using namespace std;
using namespace PLMD;

ActionWithArguments::ActionWithArguments(const ActionOptions&ao):
  Action(ao)
{
  string args;
  parse("ARG",args);

  vector<string> c=Tools::getWords(args,",");

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
            arguments.push_back(all[j]->getValue(k));
          }
        };
      } else {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(a);
        assert(action);
        if(name=="*"){
          vector<string> s=action->getValueNames();
          for(unsigned j=0;j<s.size();j++)arguments.push_back(action->getValue(s[j]));
        } else {
          assert(action->hasNamedValue(name));
          arguments.push_back(action->getValue(name));
        }
      }
    } else if(c[i]=="*"){
      std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
      for(unsigned j=0;j<all.size();j++){
        assert(all[j]->hasNamedValue(""));
        arguments.push_back(all[j]->getValue(""));
      };
    } else{
      ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(c[i]);
      assert(action);
      assert(action->hasNamedValue(""));
      arguments.push_back(action->getValue(""));
    }
  }

  if(c.size()>0){
    log.printf("  with arguments");
    for(unsigned i=0;i<c.size();i++) log.printf(" %s",c[i].c_str());
    log.printf("\n");
  }

  clearDependencies();
  for(unsigned i=0;i<arguments.size();i++) addDependency(&arguments[i]->getAction());
}




