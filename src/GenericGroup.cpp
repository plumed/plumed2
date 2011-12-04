#include "ActionRegister.h"
#include "ActionAtomistic.h"

using namespace std;

namespace PLMD {

//+PLUMEDOC GENERIC GROUP
/**
Define a group of atoms

\par Example
The following contains a static group containing atoms 1-20.  Wherever the label
of the group appears after the GROUP keyword the specified list of atom will be used
to calculate the colvar.  
\verbatim
STATIC_GROUP LABEL=label ATOMS=1-20
\endverbatim

*/
//+ENDPLUMEDOC

class GenericGroup : public ActionAtomistic {
public:
  GenericGroup(const ActionOptions&ao);
  void calculate(){};
  void apply(){};
};

PLUMED_REGISTER_ACTION(GenericGroup,"GROUP")

GenericGroup::GenericGroup(const ActionOptions&ao):
ActionAtomistic(ao)
{
  forbidKeyword("STRIDE");
  registerKeyword( 2, "ATOMS", "the numerical indexes for the set of atoms in the group");
  setNeighbourListStyle("none");
  readActionAtomistic();
  parseAtomList("ATOMS", 0 ); printAllAtoms("contains atoms");
}

} 
