#include "ActionWithVirtualAtom.h"
#include "Atoms.h"
#include "PlumedMain.h"

using namespace std;

namespace PLMD{

ActionWithVirtualAtom::ActionWithVirtualAtom(const ActionOptions&ao):
  ActionAtomistic(ao)
{
  forbidKeyword("STRIDE");
  registerKeyword( 2, "ATOMS", "the numerical indexes for the atoms involved");
}

ActionWithVirtualAtom::~ActionWithVirtualAtom(){
  plumed.getAtoms().removeVirtualAtom(this);
}

void ActionWithVirtualAtom::readActionWithVirtualAtom(){

  // Read in atom lists and so on
  readActionAtomistic();
  parseAtomList("ATOMS", 0 ); printAllAtoms("contains atoms");
  unsigned natoms=getNumberOfAtoms();

  // Setup the values
  std::vector<double> domain(2,0.0);
  readActionWithExternalArguments( 3*natoms + 9, domain );

  // Create the virtual atom
  index=plumed.getAtoms().addVirtualAtom(this);
  AtomNumber a=AtomNumber::index(index);
  log.printf("  serial associated to this virtual atom is %d\n",a.serial()); 

  // Create places to hold the values inside ActionWithValue
  addValue("x", false, true ); addValue("y", false, true ); addValue("z", false, true );  
  derivatives.resize(natoms); f.resize(natoms);
}

void ActionWithVirtualAtom::apply(){
  const Vector & forces(plumed.getAtoms().forces[index]);
  for(unsigned i=0;i<getNumberOfAtoms();++i){
     f[i]=matmul(derivatives[i],forces);
  }
  Tensor v; v.clear(); applyForces( f, v );	
}

}
