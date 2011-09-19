#include "ActionRegister.h"
#include "GenericGroup.h"

using namespace std;

namespace PLMD {

//+PLUMEDOC GENERIC GROUP
/**
Define a group of atoms

\par Syntax
\verbatim
GROUP LABEL=label ATOMS=x,y,z,...
\endverbatim
The label is associated to a group of atoms which is then automatically
expanded when used in multi-atoms options

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(GenericGroup,"GROUP")

GenericGroup::GenericGroup(const ActionOptions&ao):
  ActionAtomistic(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  vector<unsigned> a(atoms.size());
//  for(unsigned i=0;i<atoms.size();i++) a[i]=atoms[i].index();
//  plumed.getAtoms().insertGroup(getLabel(),a);
  log.printf("  of atoms ");
  for(unsigned i=0;i<atoms.size();i++) log.printf(" %d",atoms[i].serial());
  log.printf("\n");
}

GenericGroup::~GenericGroup(){
  plumed.getAtoms().removeGroup(getLabel());
}

}
