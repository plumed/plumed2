#include "Colvar.h"
#include "PlumedMain.h"

using namespace std;
using namespace PLMD;

Colvar::Colvar(const ActionOptions&ao) :
ActionAtomistic(ao)
{
  forbidKeyword("STRIDE");
}

void Colvar::readActionColvar(){
  // Check that neighbour lists have been setup properly
  checkNeighbourLists();
  // Resize stuff for applying forces
  f.resize( getNumberOfAtoms() ); forces.resize( 3*getNumberOfAtoms()+9 );
}

void Colvar::apply(){

  const unsigned nat=f.size();
  for(unsigned i=0;i<nat;i++){
    f[i][0]=0.0; f[i][1]=0.0; f[i][2]=0.0;
  }

  Tensor v; v.clear();
  for(int i=0;i<getNumberOfValues();++i){
    if( getForces( i, forces ) ){
       for(unsigned j=0;j<nat;++j){
          f[j][0]+=forces[3*j+0];
          f[j][1]+=forces[3*j+1];
          f[j][2]+=forces[3*j+2];
       }
       v(0,0)+=forces[3*nat+0];
       v(0,1)+=forces[3*nat+1];
       v(0,2)+=forces[3*nat+2];
       v(1,0)+=forces[3*nat+3];
       v(1,1)+=forces[3*nat+4];
       v(1,2)+=forces[3*nat+5];
       v(2,0)+=forces[3*nat+6];
       v(2,1)+=forces[3*nat+7];
       v(2,2)+=forces[3*nat+8];
    }
  }
  applyForces( f, v );
}
