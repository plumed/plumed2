#include "Colvar.h"
#include "PlumedMain.h"
#include <vector>
#include <string>
#include <cassert>

using namespace std;
using namespace PLMD;

Colvar::Colvar(const ActionOptions&ao) :
ActionAtomistic(ao)
{
}


//void Colvar::requestAtoms(const vector<AtomNumber> & a){
//  ActionAtomistic::requestAtoms(a);
//  setNumberOfParameters(3*a.size()+9);
//}

void Colvar::apply(){
  vector<Vector>&   f(modifyForces());
  Tensor&           v(modifyVirial());

  for(unsigned i=0;i<f.size();i++){
    f[i][0]=0.0;
    f[i][1]=0.0;
    f[i][2]=0.0;
  }
  v.clear();

  for(int i=0;i<getNumberOfValues();++i){
    if(!getValue(i)->checkForced())continue;
    const vector<double> & derivatives(getValue(i)->getDerivatives());
    const unsigned nat=f.size();
    const double force=getValue(i)->getForce();
    for(unsigned j=0;j<nat;++j){
      f[j][0]+=force*derivatives[3*j+0];
      f[j][1]+=force*derivatives[3*j+1];
      f[j][2]+=force*derivatives[3*j+2];
    }
    v(0,0)+=force*derivatives[3*nat+0];
    v(0,1)+=force*derivatives[3*nat+1];
    v(0,2)+=force*derivatives[3*nat+2];
    v(1,0)+=force*derivatives[3*nat+3];
    v(1,1)+=force*derivatives[3*nat+4];
    v(1,2)+=force*derivatives[3*nat+5];
    v(2,0)+=force*derivatives[3*nat+6];
    v(2,1)+=force*derivatives[3*nat+7];
    v(2,2)+=force*derivatives[3*nat+8];
  }
  applyForces();
}





