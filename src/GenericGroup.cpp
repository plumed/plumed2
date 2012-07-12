#ifndef __PLUMED_GenericGroup_h
#define __PLUMED_GenericGroup_h

#include "ActionRegister.h"
#include "ActionAtomistic.h"
#include "Atoms.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC GROUP
/*
Define a group of atoms so that a particular list of atoms can be referenced with a single label
in definitions of CVs or virtual atoms.

\par Examples
This command creates a group of atoms containing atoms 1,2 and 3 and assigns the label
g1 to this group.  When the label g1 appears in atom lists it is automatically expanded 
and replaced by the list of atoms (i.e. atoms 1,2 and 3). 
\verbatim
GROUP ATOMS=1,2,3 LABEL=g1
\endverbatim

*/
//+ENDPLUMEDOC

class GenericGroup:
  public ActionAtomistic
{

public:
  GenericGroup(const ActionOptions&ao);
  ~GenericGroup();
  static void registerKeywords( Keywords& keys );
  void calculate(){};
  void apply(){};
};

PLUMED_REGISTER_ACTION(GenericGroup,"GROUP")

GenericGroup::GenericGroup(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  this->atoms.insertGroup(getLabel(),atoms);
  log.printf("  of atoms ");
  for(unsigned i=0;i<atoms.size();i++) log.printf(" %d",atoms[i].serial());
  log.printf("\n");
}

void GenericGroup::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("atoms", "ATOMS", "the numerical indexes for the set of atoms in the group");
}

GenericGroup::~GenericGroup(){
  atoms.removeGroup(getLabel());
}

}

#endif
