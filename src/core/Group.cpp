/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

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
#include "Group.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "Value.h"
#include "ActionWithValue.h"
#include "ActionWithVirtualAtom.h"
#include "tools/IFile.h"
#include "tools/Tools.h"
#include <string>
#include <vector>
#include <algorithm>

namespace PLMD {

//+PLUMEDOC GENERIC GROUP
/*
Define a group of atoms so that a particular list of atoms can be referenced with a single label in definitions of CVs or virtual atoms.

The GROUP command can be used to define a list of atoms so that the group can be referenced in the definitions of other
CVs or virtual atoms as shown below:

```plumed
o: GROUP ATOMS=1,4,7,11,14
h: GROUP ATOMS=2,3,5,6,8,9,12,13
# compute the coordination among the two groups
c: COORDINATION GROUPA=o GROUPB=h R_0=0.3
# same could have been obtained without GROUP, just writing:
# c: COORDINATION GROUPA=1,4,7,11,14 GROUPB=2,3,5,6,8,9,12,13

# print the coordination on file 'colvar'
PRINT ARG=c FILE=colvar
```

The first group command here creates a group of atoms called `o` containing atoms 1, 4, 7, 11 and 14.
The second group command here creates a group of toms called `h` containing atoms 2, 3, 5, 6, 8, 9, 12, and 13.

As discussed on [this page](specifying_atoms.md) atoms in groups can be listed as comma separated numbers (i.e. `1,2,3,10,45,7,9`),
simple positive ranges (i.e. `20-40`), ranges with a stride either positive or negative (i.e. `20-40:2` or `80-50:-2`) or as
comma separated combinations of all the former methods (`1,2,4,5,10-20,21-40:2,80-50:-2`).

Oftentime people will store the definitions in a separate file as has been done in the following example input.

```plumed
INCLUDE FILE=extras/groups.dat
# compute the coordination among the two groups
c: COORDINATION GROUPA=groupa GROUPB=groupb R_0=0.3
# print the coordination on file 'colvar'
PRINT ARG=c FILE=colvar
```

Storing the groups in the extra file is particularly useful if the groups include list of thousand atoms. Putting these definitions
in a separate file ensures that the main plumed.dat file is not cluttered.  You can even use a [GROMACS index file](https://manual.gromacs.org/archive/5.0.4/online/ndx.html)
to hold the groups as illustrated in the input below:

```plumed
# import group named 'Protein' from file index.ndx
pro: GROUP NDX_FILE=extras/index.ndx NDX_GROUP=Protein
# dump all the atoms of the protein on a trajectory file
DUMPATOMS ATOMS=pro FILE=traj.gro
```

Notice that you use the keyword `NDX_FILE` to set the name of the index file and `NDX_GROUP` to set the name of the group to be imported (default is first one).
Further notice that starting from version 2.10 it is possible to directly use an `@ndx:` selector in the input to an action as shwn below:

```plumed
DUMPATOMS ATOMS={@ndx:{extras/index.ndx Protein}} FILE=traj.gro
```

Notice that it is possible to remove atoms from the list of atoms specified in a GROUP by using the keyword `REMOVE` as shown below:

```plumed
# take one atom every three, that is oxygens
ox: GROUP ATOMS=1-90:3
# take the remaining atoms, that is hydrogens
hy: GROUP ATOMS=1-90 REMOVE=ox
DUMPATOMS ATOMS=ox FILE=ox.gro
DUMPATOMS ATOMS=hy FILE=hy.gro
```

You can also `SORT` the atoms in a group into ascending order by using the SORT keyword as shown below:

```plumed
# If you ask for this group in the input to another action the atoms will be
# in ascending order i.e. the specified atoms will be 2,3,4,4,5,6
g: GROUP ATOMS=5,4,6,3,4,2 SORT
```

or you can sort the atoms and remove duplicated atoms by using the `UNIQUE` flag

```plumed
# If you ask for this group in the input to another action the duplicate atom specifications will be removed
# and the atoms will be ascending order i.e. the specified atoms will be 2,3,4,5,6
g: GROUP ATOMS=5,4,6,3,4,2 UNIQUE
```

When you used the GROUP command the flow as follows:

- If `ATOMS` is present, then take the ordered list of atoms from the `ATOMS` keyword as a starting list.
- Alternatively, if `NDX_FILE` is present, use the list obtained from the gromacs group.
- If `REMOVE` is present, then remove the first occurrence of each of these atoms from the list.
  If one tries to remove an atom that was not listed plumed adds a notice in the output.
  An atom that is present twice in the original list might be removed twice.
- If `SORT` is present, then the resulting list is sorted by increasing serial number.
- If `UNIQUE` is present, then the resulting list is sorted by increasing serial number _and_ duplicate elements are removed.

Notice that this command just creates a shortcut, and does not imply any real calculation.
So, having a huge group defined does not slow down your calculation in any way.
It is just convenient to better organize input files.

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(Group,"GROUP")

Group::Group(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao) {
  parseAtomList("ATOMS",atoms);
  std::string ndxfile,ndxgroup;
  parse("NDX_FILE",ndxfile);
  parse("NDX_GROUP",ndxgroup);
  if(ndxfile.length()>0 && atoms.size()>0) {
    error("either use explicit atom list or import from index file");
  }
  if(ndxfile.length()==0 && ndxgroup.size()>0) {
    error("NDX_GROUP can be only used is NDX_FILE is also used");
  }

  if(ndxfile.length()>0) {

    std::vector<AtomNumber> add;
    std::vector<std::string> words;
    words.emplace_back("@ndx: " + ndxfile + " " + ndxgroup);
    interpretAtomList(words,add);
    atoms.insert(atoms.end(),add.begin(),add.end());
  }

  std::vector<AtomNumber> remove;
  parseAtomList("REMOVE",remove);
  if(remove.size()>0) {
    std::vector<AtomNumber> notfound;
    unsigned k=0;
    log<<"  removing these atoms from the list:";
    for(unsigned i=0; i<remove.size(); i++) {
      const auto it = find(atoms.begin(),atoms.end(),remove[i]);
      if(it!=atoms.end()) {
        if(k%25==0) {
          log<<"\n";
        }
        log<<" "<<(*it).serial();
        k++;
        atoms.erase(it);
      } else {
        notfound.push_back(remove[i]);
      }
    }
    log<<"\n";
    if(notfound.size()>0) {
      log<<"  the following atoms were not found:";
      for(unsigned i=0; i<notfound.size(); i++) {
        log<<" "<<notfound[i].serial();
      }
      log<<"\n";
    }
  }

  bool sortme=false;
  parseFlag("SORT",sortme);
  if(sortme) {
    log<<"  atoms are sorted\n";
    sort(atoms.begin(),atoms.end());
  }
  bool uniqueFlag=false;
  parseFlag("UNIQUE",uniqueFlag);
  if(uniqueFlag) {
    log<<"  sorting atoms and removing duplicates\n";
    Tools::removeDuplicates(atoms);
  }

  log.printf("  list of atoms:");
  for(unsigned i=0; i<atoms.size(); i++) {
    if(i%25==0) {
      log<<"\n";
    }
    log<<" "<<atoms[i].serial();
  }
  log.printf("\n");
}

void Group::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("atoms", "ATOMS", "the numerical indexes for the set of atoms in the group");
  keys.add("atoms", "REMOVE","remove these atoms from the list");
  keys.addFlag("SORT",false,"sort the resulting list");
  keys.addFlag("UNIQUE",false,"sort atoms and remove duplicated ones");
  keys.add("optional", "NDX_FILE", "the name of index file (gromacs syntax)");
  keys.add("optional", "NDX_GROUP", "the name of the group to be imported (gromacs syntax) - first group found is used by default");
}

std::vector<std::string> Group::getGroupAtoms() const {
  std::vector<std::string> atoms_str(atoms.size());
  for(unsigned i=0; i<atoms.size(); ++i) {
    std::pair<std::size_t,std::size_t> a = getValueIndices( atoms[i] );
    if( xpos[a.first]->getNumberOfValues()==1 ) {
      atoms_str[i] = (xpos[a.first]->getPntrToAction())->getLabel();
    } else {
      Tools::convert( atoms[i].serial(), atoms_str[i] );
    }
  }
  return atoms_str;
}

}

