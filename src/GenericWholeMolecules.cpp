#include "ActionAtomistic.h"
#include "ActionPilot.h"
#include "ActionRegister.h"
#include "Vector.h"
#include "AtomNumber.h"
#include "Tools.h"

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


class GenericWholeMolecules:
  public ActionPilot,
  public ActionAtomistic
{
  vector<vector<AtomNumber> > groups;
  Vector & modifyPosition(AtomNumber);
public:
  GenericWholeMolecules(const ActionOptions&ao);
  void calculate();
  void apply(){};
};

PLUMED_REGISTER_ACTION(GenericWholeMolecules,"WHOLEMOLECULES")

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
    string is; Tools::convert(i,is);
    string name="GROUP"+is;
    vector<AtomNumber> group;
    parseAtomList(name,group);
    if(group.size()==0)break;
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

