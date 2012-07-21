#include "ActionSet.h"

using namespace std;
using namespace PLMD;

ActionSet::ActionSet(PlumedMain&p):
plumed(p){
}

ActionSet::~ActionSet()
{
  for(int i=size()-1;i>=0;i--) delete (*this)[i];
}

void ActionSet::clearDelete(){
  for(int i=size()-1;i>=0;i--) delete (*this)[i];
  clear();
}


std::string ActionSet::getLabelList() const{
  std::string outlist;
  for(const_iterator p=begin();p!=end();++p){
    outlist+=dynamic_cast<Action*>(*p)->getLabel()+" ";
  };
  return  outlist;
}


