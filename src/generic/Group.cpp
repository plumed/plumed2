/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

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

#include "core/ActionRegister.h"
#include "core/ActionAtomistic.h"
#include "core/Atoms.h"

using namespace std;

namespace PLMD{
namespace generic{

//+PLUMEDOC GENERIC GROUP
/*
Define a group of atoms so that a particular list of atoms can be referenced with a single label
in definitions of CVs or virtual atoms.

Notice that this command just creates a shortcut, and does not imply any real calculation.
It is just convenient to better organize input files. Might be used in combination with
the \ref INCLUDE command so as to store long group definitions in a separate files.

\par Examples

This command create a group of atoms containing atoms 1,4,7,11 and 14 (labeled 'o'), and another containing
atoms 2,3,5,6,8,9,12,13 (labeled 'h'):
\verbatim
o: GROUP ATOMS=1,4,7,11,14
h: GROUP ATOMS=2,3,5,6,8,9,12,13
# compute the coordination among the two groups
c: COORDINATION GROUPA=o GROUPB=h R_0=0.3

# same could have been obtained without GROUP, just writing:
# c: COORDINATION GROUPA=1,4,7,11,14 GROUPB=2,3,5,6,8,9,12,13
\endverbatim
(see also \ref COORDINATION)

Groups can be conveniently stored in a separate file
\verbatim
INCLUDE groups.dat
c: COORDINATION GROUPA=o GROUPB=h R_0=0.3
\endverbatim
(see also \ref INCLUDE and \ref COORDINATION).
The groups.dat file could be very long and include lists of thousand atoms without cluttering the main plumed.dat file.

*/
//+ENDPLUMEDOC

class Group:
  public ActionAtomistic
{

public:
  Group(const ActionOptions&ao);
  ~Group();
  static void registerKeywords( Keywords& keys );
  void calculate(){}
  void apply(){}
};

PLUMED_REGISTER_ACTION(Group,"GROUP")

Group::Group(const ActionOptions&ao):
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

void Group::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("atoms", "ATOMS", "the numerical indexes for the set of atoms in the group");
}

Group::~Group(){
  atoms.removeGroup(getLabel());
}

}
}

