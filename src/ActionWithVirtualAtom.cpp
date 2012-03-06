#include "ActionWithVirtualAtom.h"
#include "Atoms.h"
#include "PlumedMain.h"

using namespace std;

namespace PLMD{

void ActionWithVirtualAtom::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  keys.add("atoms","ATOMS","the list of atoms which are involved the virtual atom's definition");
  keys.addFlag("CALC_GRADIENTS", false, "  calculate the vector of gradients");
}

ActionWithVirtualAtom::ActionWithVirtualAtom(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao)
{
  index=atoms.addVirtualAtom(this);
  log.printf("  serial associated to this virtual atom is %d\n",index.serial());
}

ActionWithVirtualAtom::~ActionWithVirtualAtom(){
  atoms.removeVirtualAtom(this);
}

void ActionWithVirtualAtom::apply(){
  const Vector & f(atoms.forces[index.index()]);
  for(unsigned i=0;i<getNumberOfAtoms();i++) modifyForces()[i]=matmul(derivatives[i],f);
}

void ActionWithVirtualAtom::requestAtoms(const std::vector<AtomNumber> & a){
  ActionAtomistic::requestAtoms(a);
  derivatives.resize(a.size());
}

void ActionWithVirtualAtom::setGradients(){
  gradients.clear();
  for(unsigned i=0;i<getNumberOfAtoms();i++){
    AtomNumber an=getAbsoluteIndex(i);
    // this case if the atom is a virtual one 	 
    if(atoms.isVirtualAtom(an)){
      const ActionWithVirtualAtom* a=atoms.getVirtualAtomsAction(an);
      for(std::map<AtomNumber,Tensor>::const_iterator p=a->gradients.begin();p!=a->gradients.end();++p){
        gradients[(*p).first]+=matmul(derivatives[i],(*p).second);
      }
    // this case if the atom is a normal one 	 
    } else {
      gradients[an]+=derivatives[i];
    }
  }
}


void ActionWithVirtualAtom::setGradientsIfNeeded(){
  if(isOptionOn("GRADIENTS")) { 
	setGradients() ;	
  }
}

}
