#include "ActionSetup.h"
#include "ActionRegister.h"
#include "Vector.h"
#include "AtomNumber.h"
#include "Tools.h"
#include "Atoms.h"
#include "PlumedMain.h"

#include <vector>
#include <string>

using namespace std;
using namespace PLMD;

namespace PLMD {

//+PLUMEDOC GENERIC WHOLEMOLECULES
/**
Rebuild molecules with pbc

\par syntax
\verbatim
WHOLEMOLECULES [STRIDE=s] GROUP0=list0 [ GROUP1=list1 [ GROUP2=list2 [ ... ] ] ]
\endverbatim

Similar to the ALIGN_ATOMS keyword of plumed 1. It rebuilds molecules
correctly according to pbc. It can rebuild multiple groups, and
it can be applied more than once. To see its effects, use
a variable without pbc or the \ref DUMPATOMS directive.

\attention
This directive is modifying the stored position in the precise moment
when it get executed. This means that only collective variables
which are below this in the input script will see the corrected positions.
As a general rule, put it at the top of the input file. Also, unless you
know exactly what you are doing, leave the default stride (1), so that
it acts at every step.
*/
//+ENDPLUMEDOC


class GenericWholeMolecules : public ActionSetup {
public:
  GenericWholeMolecules(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(GenericWholeMolecules,"WHOLEMOLECULES")

GenericWholeMolecules::GenericWholeMolecules(const ActionOptions&ao):
ActionSetup(ao)
{
  registerKeyword(2,"MOLECULE","the atoms that make up a molecule that you wish to align. To specify multiple molecules use a list of MOLECULE keywords: MOLECULE1, MOLECULE2,...");
  allowKeyword("MOLECULE"); 
  readActionSetup();

  if ( !testForKey("MOLECULE") ) error("no molecules specified in input to WHOLEMOLECULES");

  std::vector<std::string> strings; Atoms& atoms(plumed.getAtoms());
  std::string myname;
  if( testForNumberedKeys("MOLECULE") ){
    std::string num; 
    for(int i=1;;++i ){
       Tools::convert(i,num); myname="chain" + num;
       if( !parseNumberedVector( "MOLECULE", i, strings ) ) break;
       atoms.addMolecule( *this, myname, strings );
       log.printf("  molecule %d contains the following atoms : ",i);
       for(unsigned i=0;i<strings.size();++i) log.printf("%s ", strings[i].c_str() );
       log.printf("\n");
     }
  } else {
     parseVector("MOLECULE",strings); myname="chain1";
     atoms.addMolecule( *this, myname, strings );
     log.printf("  molecule 1 contains the following atoms : ");
     for(unsigned i=0;i<strings.size();++i) log.printf("%s ", strings[i].c_str() );
     log.printf("\n");
  }  
}

}

