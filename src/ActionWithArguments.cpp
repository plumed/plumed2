/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ActionWithArguments.h"
#include "ActionWithValue.h"
#include "PlumedMain.h"
#include "ActionSet.h"

using namespace std;
using namespace PLMD;

void ActionWithArguments::registerKeywords(Keywords& keys){
  keys.reserve("compulsory","ARG","the input for this action is the output from one or more other actions. The particular output that you used is referenced using that action of interests label. If the label appears on its own then the value of the relevant Action is taken.  If * or *.* appears the information from all arguments is taken.  Some actions have multi-component outputs, each component of the output has a specific label so for instance an action labelled dist may have three componets x, y and z.  To take just the x component you should use dist.x, if you wish to take all three components then use dist.*");
}

void ActionWithArguments::parseArgumentList(const std::string&key,std::vector<Value*>&arg){
  vector<string> c; arg.clear(); parseVector(key,c);

  for(unsigned i=0;i<c.size();i++){
      std::size_t dot=c[i].find_first_of('.');
      string a=c[i].substr(0,dot);
      string name=c[i].substr(dot+1);
      if(c[i].find(".")!=string::npos){    // if it contains a dot:
        if(a=="*" && name=="*"){
           // Take all values from all actions
           std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
           if( all.empty() ) error("your input file is not telling plumed to calculate anything");
           for(unsigned j=0;j<all.size();j++){
             for(int k=0;k<all[j]->getNumberOfComponents();++k) arg.push_back(all[j]->copyOutput(k));
           }
        } else if ( name=="*"){
           // Take all the values from an action with a specific name
           ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(a);
           if(!action) error("cannot find action named " + a);
           for(int k=0;k<action->getNumberOfComponents();++k) arg.push_back(action->copyOutput(k));
        } else if ( a=="*" ){
           // Take components from all actions with a specific name
           std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
           if( all.empty() ) error("your input file is not telling plumed to calculate anything");
           std::string lab; unsigned nval=0;
           for(unsigned j=0;j<all.size();j++){
              std::string flab; flab=all[j]->getLabel() + "." + name;
              if( all[j]->exists(flab) ){ arg.push_back(all[j]->copyOutput(flab)); nval++; }
           }
           if(nval==0) error("found no actions with a component called " + name );
        } else {
           // Take values with a specific name
           ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(a);
           if(!action) error("cannot find action named " + a);
           if( !(action->exists(c[i])) ) error("action " + a + " has no component named " + name );
           arg.push_back(action->copyOutput(c[i]));
        }
      } else {    // if it doesn't contain a dot
        if(c[i]=="*"){
           // Take all values from all actions
           std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
           if( all.empty() ) error("your input file is not telling plumed to calculate anything");
           for(unsigned j=0;j<all.size();j++){
             for(int k=0;k<all[j]->getNumberOfComponents();++k) arg.push_back(all[j]->copyOutput(k));
           }
        } else {
           ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(c[i]);
           if(!action) error("cannot find action named " + c[i]);
           if( !(action->exists(c[i])) ) error("action " + c[i] + " has no component named " + c[i] );
           arg.push_back(action->copyOutput(c[i]));
        }
      }
  }
}

void ActionWithArguments::requestArguments(const vector<Value*> &arg){
  plumed_massert(!lockRequestArguments,"requested argument list can only be changed in the prepare() method");
  arguments=arg;
  clearDependencies();
  std::string fullname,name;
  for(unsigned i=0;i<arguments.size();i++){
     fullname=arguments[i]->getName();
     if(fullname.find(".")!=string::npos){
       std::size_t dot=fullname.find_first_of('.');
       name=fullname.substr(0,dot);
     } else {
       name=fullname;
     }
     ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>(name);
     plumed_massert(action,"cannot find action named (in requestArguments - this is weird)" + name);
     addDependency(action);
  }
}

ActionWithArguments::ActionWithArguments(const ActionOptions&ao):
  Action(ao),
  lockRequestArguments(false)
{
  if( keywords.exists("ARG") ){
     vector<Value*> arg;
     parseArgumentList("ARG",arg);

     if(!arg.empty()){
       log.printf("  with arguments");
       for(unsigned i=0;i<arg.size();i++) log.printf(" %s",arg[i]->getName().c_str());
       log.printf("\n");
     }
     requestArguments(arg);
  }
}

void ActionWithArguments::calculateNumericalDerivatives( ActionWithValue* a ){
  if(!a){
    a=dynamic_cast<ActionWithValue*>(this);
    plumed_massert(a,"cannot compute numerical derivatives for an action without values");
  }

  const int nval=a->getNumberOfComponents();
  const int npar=arguments.size();
  std::vector<double> value (nval*npar);
  for(int i=0;i<npar;i++){
    double arg0=arguments[i]->get();
    arguments[i]->set(arg0+sqrt(epsilon));
    a->calculate();
    arguments[i]->set(arg0);
    for(unsigned j=0;j<nval;j++){
      value[i*nval+j]=a->getOutputQuantity(j);
    }
  }
  a->calculate();
  a->clearDerivatives();
  std::vector<double> value0(nval);
  for(unsigned j=0;j<nval;j++){
    Value* v=a->copyOutput(j);
    if( v->hasDerivatives() ) for(int i=0;i<npar;i++) v->addDerivative(i,(value[i*nval+j]-a->getOutputQuantity(j))/sqrt(epsilon));
  }
}

double ActionWithArguments::getProjection(unsigned i,unsigned j)const{
  plumed_massert(i<arguments.size()," making projections with an index which  is too large");
  plumed_massert(j<arguments.size()," making projections with an index which  is too large");
  const Value* v1=arguments[i];
  const Value* v2=arguments[j];
  return Value::projection(*v1,*v2);
}







