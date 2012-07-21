/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
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
