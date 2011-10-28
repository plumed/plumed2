#include "ActionSetup.h"
#include "ActionRegister.h"
#include "Vector.h"
#include "AtomNumber.h"
#include "Tools.h"
#include "Atoms.h"
#include "PlumedMain.h"

using namespace std;
using namespace PLMD;

namespace PLMD {

//+PLUMEDOC TOPOLOGY PROTEIN_TOPOLOGY
/**
This command calls on plumed to read a pdb file.  This is then stored and features
like the atoms in the backbone or in specific residues can be used as short cuts
to define the atoms that are involved in the various colvars. 

\par Example

*/
//+ENDPLUMEDOC

class TopologyProtein : public ActionSetup {
friend class Atoms;
public:
  TopologyProtein(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(TopologyProtein,"PROTEIN_TOPOLOGY")

TopologyProtein::TopologyProtein(const ActionOptions&ao):
ActionSetup(ao)
{
  registerKeyword(1,"FILE","a pdb file containing the protein topology");
  readActionSetup(); 

  std::string filen; parse("FILE",filen);  
  log.printf("  reading topology from file %s\n",filen.c_str() );
  plumed.getAtoms().readTopology( *this, "protein", filen );
}

}

