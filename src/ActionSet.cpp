#include "ActionSet.h"
#include "PlumedMain.h"

using namespace std;
using namespace PLMD;

ActionSet::ActionSet(PlumedMain&p):
plumed(p){
}

void ActionSet::clearDelete(){
  for(int i=size()-1;i>=0;i--) delete (*this)[i];
  clear();
}

