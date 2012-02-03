#include "ActionAtomistic.h"
#include "ActionPilot.h"
#include "ActionRegister.h"
#include "Vector.h"
#include "AtomNumber.h"
#include "Tools.h"
#include "PlumedMain.h"
#include "Atoms.h"

#include <vector>
#include <string>

using namespace std;
using namespace PLMD;

namespace PLMD {

//+PLUMEDOC GENERIC WHOLEMOLECULES
/**
This action is used to rebuild molecules that can become split by the periodic
boundary conditions in a manner similar to the ALIGN_ATOMS keyword of plumed1.

Running some CVs without this command can cause there to be discontinuities changes
in the CV value and artifacts in the calculations.  This command can be applied 
more than once.  To see what effect is has use a variable without pbc or use
the \ref DUMPATOMS directive to output the atomic positions.

\attention
This directive modifies the stored position at the precise moment
it is executed. This means that only collective variables
which are below it in the input script will see the corrected positions.
As a general rule, put it at the top of the input file. Also, unless you
know exactly what you are doing, leave the default stride (1), so that
this action is performed at every MD step.

\par Examples
This command instructs plumed to reconstruct the molecule containing atoms 1-20
at every step of the calculation.

\verbatim
WHOLEMOLECULES STRIDE=1 MOLECULE=1-20
\endverbatim

*/
//+ENDPLUMEDOC


class GenericWholeMolecules:
  public ActionPilot,
  public ActionAtomistic
{
  vector<vector<AtomNumber> > groups;
  Vector & modifyPosition(AtomNumber);
public:
  GenericWholeMolecules(const ActionOptions&ao);
  static void registerKeywords( Keywords& keys );
  void calculate();
  void apply(){};
};

PLUMED_REGISTER_ACTION(GenericWholeMolecules,"WHOLEMOLECULES")

void GenericWholeMolecules::registerKeywords( Keywords& keys ){
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which molecules are reassembled.  Unless you are completely certain about what you are doing leave this set equal to 1!");
  keys.add("input","MOLECULE","the atoms that make up a molecule that you wish to align. To specify multiple molecules use a list of MOLECULE keywords: MOLECULE1, MOLECULE2,...");
}

inline
Vector & GenericWholeMolecules::modifyPosition(AtomNumber i){
  return plumed.getAtoms().positions[i.index()];
}

GenericWholeMolecules::GenericWholeMolecules(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionAtomistic(ao)
{
  vector<AtomNumber> merge;
  for(int i=0;;i++){
    //string is; Tools::convert(i,is);
    //string name="MOLECULE"+is;
    vector<AtomNumber> group;
    if(!parseNumberedAtomList("MOLECULE",i,group) ) break;
    groups.push_back(group);
    merge.insert(merge.end(),group.begin(),group.end());
  }
  checkRead();
  Tools::removeDuplicates(merge);
  requestAtoms(merge);
}

void GenericWholeMolecules::calculate(){
  for(unsigned i=0;i<groups.size();++i){
    for(unsigned j=0;j<groups[i].size()-1;++j){
      Vector & first (modifyPosition(groups[i][j]));
      Vector & second (modifyPosition(groups[i][j+1]));
      second=first+pbcDistance(first,second);
    }
  }
}



}

