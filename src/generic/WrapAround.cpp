/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Vector.h"
#include "tools/AtomNumber.h"
#include "tools/Tools.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/GenericMolInfo.h"

#include <vector>

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC WRAPAROUND
/*
Rebuild periodic boundary conditions around chosen atoms.

This action modifies the position of the atoms indicated by ATOMS by shifting them by lattice vectors so that they are
as close as possible to the atoms indicated by AROUND. More precisely, for every atom i
in the ATOMS list the following procedure is performed:
- The atom j among those in the AROUND list is searched that is closest to atom i.
- The atom i is replaced with its periodic image that is closest to atom j.

This action works similarly to [WHOLEMOLECULES](WHOLEMOLECULES.md) in that it replaces atoms coordinate. Notice that only
atoms specified with ATOMS are replaced, and that, at variance with [WHOLEMOLECULES](WHOLEMOLECULES.md),
the order in which atoms are specified is irrelevant.

This is often convenient at a post processing stage (using the driver), but sometime
it is required during the simulation if collective variables need atoms to be in a specific periodic image.

!!! caution "modifies stored positions

    This directive modifies the stored position at the precise moment it is executed. This means that only collective variables which are below it in
    the input script will see the corrected positions. As a general rule, put it at the top of the input file. Also, unless you know exactly what you are doing,
    leave the default stride (1), so that this action is performed at every MD step.

The computational cost of this action grows with the product
of the size of the two lists (ATOMS and AROUND), so this action can become very expensive.
If you are using it to analyze a trajectory this is usually not a big problem. If you use it to
analyze a simulation on the fly, e.g. with [DUMPATOMS](DUMPATOMS.md) to store a properly wrapped trajectory,
consider using the STRIDE keyword here (with great care).


## Examples

This command instructs plumed to move all the ions to their periodic image that is as close as possible to
the rna group.

```plumed
rna: GROUP ATOMS=1-100
ions: GROUP ATOMS=101-110
# first make the rna molecule whole
WHOLEMOLECULES ENTITY0=rna
WRAPAROUND ATOMS=ions AROUND=rna
DUMPATOMS FILE=dump.xyz ATOMS=rna,ions
```

In case you want to do it during a simulation and you only care about wrapping the ions in
the `dump.xyz` file, you can use the following input:

```plumed
# add some restraint that do not require molecules to be whole:
a: TORSION ATOMS=1,2,10,11
RESTRAINT ARG=a AT=0.0 KAPPA=5


# then do the things that are required for dumping the trajectory
# notice that they are all done every 100 steps, so as not to
# unnecessarily overload the calculation

rna: GROUP ATOMS=1-100
ions: GROUP ATOMS=101-110
# first make the rna molecule whole
WHOLEMOLECULES ENTITY0=rna STRIDE=100
WRAPAROUND ATOMS=ions AROUND=rna STRIDE=100
DUMPATOMS FILE=dump.xyz ATOMS=rna,ions STRIDE=100
```

Notice that if the biased variable requires a molecule to be whole, you might have to put
the [WHOLEMOLECULES](WHOLEMOLECULES.md) command before computing that variable and leave the default STRIDE=1.

This command instructs plumed to center all atoms around the center of mass of a solute molecule.

```plumed
solute: GROUP ATOMS=1-100
all: GROUP ATOMS=1-1000
# center of the solute:
# notice that since plumed 2.2 this also works if the
# solute molecule is broken
com: COM ATOMS=solute
# notice that we wrap around a single atom. this should be fast
WRAPAROUND ATOMS=all AROUND=com
DUMPATOMS FILE=dump.xyz ATOMS=all
```

Notice that whereas [WHOLEMOLECULES](WHOLEMOLECULES.md) is designed to make molecules whole,
WRAPAROUND can easily break molecules. In the last example,
if solvent (atoms 101-1000) is made e.g. of water, then water
molecules could be broken by WRAPAROUND (hydrogen could end up
in an image and oxygen in another one).
One solution is to use [WHOLEMOLECULES](WHOLEMOLECULES.md) on _all_ the water molecules
after WRAPAROUND. This is tedious. A better solution is to use the
GROUPBY option which is going
to consider the atoms listed in ATOMS as a list of groups
each of size GROUPBY. The first atom of the group will be brought
close to the AROUND atoms. The following atoms of the group
will be just brought close to the first atom of the group.
Assuming that oxygen is the first atom of each water molecules,
in the following examples all the water oxygen atoms will be brought
close to the solute, and all the hydrogen atoms will be kept close
to their related oxygen.

```plumed
solute: GROUP ATOMS=1-100
water: GROUP ATOMS=101-1000
com: COM ATOMS=solute
# notice that we wrap around a single atom. this should be fast
WRAPAROUND ATOMS=solute AROUND=com
# notice that we wrap around a single atom. this should be fast
WRAPAROUND ATOMS=water AROUND=com GROUPBY=3
DUMPATOMS FILE=dump.xyz ATOMS=solute,water
```

*/
//+ENDPLUMEDOC


class WrapAround:
  public ActionPilot,
  public ActionAtomistic {
  // cppcheck-suppress duplInheritedMember
  std::vector<Vector> refatoms;
  std::vector<std::pair<std::size_t,std::size_t> > p_atoms;
  std::vector<std::pair<std::size_t,std::size_t> > p_reference;
  unsigned groupby;
  bool pair_;
public:
  explicit WrapAround(const ActionOptions&ao);
  static void registerKeywords( Keywords& keys );
  bool actionHasForces() override {
    return false;
  }
  void calculate() override;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(WrapAround,"WRAPAROUND")

void WrapAround::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which molecules are reassembled.  Unless you are completely certain about what you are doing leave this set equal to 1!");
  keys.add("atoms","AROUND","reference atoms");
  keys.add("atoms","ATOMS","wrapped atoms");
  keys.add("compulsory","GROUPBY","1","group atoms so as not to break molecules");
  keys.addFlag("PAIR", false, "Pair atoms in AROUND and ATOMS groups");
}

WrapAround::WrapAround(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  groupby(1),
  pair_(false) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  std::vector<AtomNumber> reference;
  parseAtomList("AROUND",reference);
  parse("GROUPBY",groupby);
  parseFlag("PAIR", pair_);

  log.printf("  atoms in reference :");
  for(unsigned j=0; j<reference.size(); ++j) {
    log.printf(" %d",reference[j].serial() );
  }
  log.printf("\n");
  log.printf("  atoms to be wrapped :");
  for(unsigned j=0; j<atoms.size(); ++j) {
    log.printf(" %d",atoms[j].serial() );
  }
  log.printf("\n");
  if(groupby>1) {
    log<<"  atoms will be grouped by "<<groupby<<"\n";
  }
  if(pair_) {
    log.printf("  pairing atoms and references\n");
  }

  if(atoms.size()%groupby!=0) {
    error("number of atoms should be a multiple of groupby option");
  }
  // additional checks with PAIR
  if(pair_ && atoms.size()!=reference.size()*groupby) {
    error("with PAIR you must have: #ATOMS = #AROUND * #GROUPBY");
  }

  checkRead();

  // do not remove duplicates with pair
  if(!pair_) {
    if(groupby<=1) {
      Tools::removeDuplicates(atoms);
    }
    Tools::removeDuplicates(reference);
  }

  std::vector<AtomNumber> merged(atoms.size()+reference.size());
  merge(atoms.begin(),atoms.end(),reference.begin(),reference.end(),merged.begin());
  p_atoms.resize( atoms.size() );
  for(unsigned i=0; i<atoms.size(); ++i) {
    p_atoms[i] = getValueIndices( atoms[i] );
  }
  refatoms.resize( reference.size() );
  p_reference.resize( reference.size() );
  for(unsigned i=0; i<reference.size(); ++i) {
    p_reference[i] = getValueIndices( reference[i] );
  }
  Tools::removeDuplicates(merged);
  requestAtoms(merged);
  doNotRetrieve();
  doNotForce();
}

void WrapAround::calculate() {
  for(unsigned j=0; j<p_reference.size(); ++j) {
    refatoms[j] = getGlobalPosition(p_reference[j]);
  }

  for(unsigned i=0; i<p_atoms.size(); i+=groupby) {
    Vector second, first=getGlobalPosition(p_atoms[i]);
    double mindist2=std::numeric_limits<double>::max();
    int closest=-1;
    if(pair_) {
      closest = i/groupby;
    } else {
      for(unsigned j=0; j<p_reference.size(); ++j) {
        second=refatoms[j];
        const Vector distance=pbcDistance(first,second);
        const double distance2=modulo2(distance);
        if(distance2<mindist2) {
          mindist2=distance2;
          closest=j;
        }
      }
      plumed_massert(closest>=0,"closest not found");
    }
    second=refatoms[closest];
// place first atom of the group
    first=second+pbcDistance(second,first);
    setGlobalPosition(p_atoms[i],first);
// then place other atoms close to the first of the group
    for(unsigned j=1; j<groupby; j++) {
      second=getGlobalPosition(p_atoms[i+j]);
      setGlobalPosition( p_atoms[i+j], first+pbcDistance(first,second) );
    }
  }
}



}

}
