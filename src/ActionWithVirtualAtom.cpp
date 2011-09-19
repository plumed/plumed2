#include "ActionWithVirtualAtom.h"
#include "Atoms.h"
#include "PlumedMain.h"

using namespace std;

namespace PLMD{

ActionWithVirtualAtom::ActionWithVirtualAtom(const ActionOptions&ao):
  ActionAtomistic(ao)
{
  index=plumed.getAtoms().addVirtualAtom(this);
  AtomNumber a=AtomNumber::index(index);
  log.printf("  serial associated to this virtual atom is %d\n",a.serial());
}

ActionWithVirtualAtom::~ActionWithVirtualAtom(){
  plumed.getAtoms().removeVirtualAtom(this);
}

void ActionWithVirtualAtom::apply(){
  const Vector & f(plumed.getAtoms().forces[index]);
  for(unsigned i=0;i<getNatoms();i++) modifyForces()[i]=matmul(derivatives[i],f);
  applyForces();	
}

//void ActionWithVirtualAtom::requestAtoms(const std::vector<AtomNumber> & a){
//  ActionAtomistic::requestAtoms(a);
//  derivatives.resize(a.size());
//}

}
