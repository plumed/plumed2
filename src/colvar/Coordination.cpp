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
#include "CoordinationBase.h"
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR COORDINATION
/*
Calculate coordination numbers.

This keyword can be used to calculate the number of contacts between two groups of atoms
and is defined as

$$
\sum_{i\in A} \sum_{i\in B} s_{ij}
$$

where $s_{ij}$ is 1 if the contact between atoms $i$ and $j$ is formed and zero otherwise.
In actuality, $s_{ij}$ is replaced with a switching function so as to ensure that the calculated CV has continuous derivatives.

The default switching function is:

$$
s_{ij} = \frac{ 1 - \left(\frac{r_{ij}-d_0}{r_0}\right)^n } { 1 - \left(\frac{r_{ij}-d_0}{r_0}\right)^m }
$$

but it can be changed using the optional SWITCH option.  You can find more information about the various switching functions that you can
use with this action in the documentation for [LESS_THAN](LESS_THAN.md).

To make your __calculation faster you can use a neighbor list__: every NL_STRIDE steps the pairs closer than NL_CUTOFF are listed,
and only those are computed at the following steps.

- With NLIST the list is built by checking all the pairs. Use it with small groups (up to about a hundred atoms), or when one of the two groups has only a few atoms (fewer than about ten).
- With NLISTCELLS the list is built with link cells, at a cost that grows with the number of atoms rather than with the number of pairs. Use it for larger groups:
  building the list is about 4 times faster with 1000 atoms, and about 30 times with 10000. The pairs are the same, in another order.
  When one of the two groups has at most 8 atoms, NLISTCELLS builds the list as NLIST: sorting all the atoms into cells would cost more than one distance per pair.
- NL_CUTOFF must be larger than D_MAX of the switching function (set D_MAX: without it the switching function is not zero at any distance, and the pairs beyond NL_CUTOFF are ignored).
  The difference NL_CUTOFF-D_MAX is a skin: the list is exact as long as no pair of atoms gets closer by more than the skin between two updates.
  For atoms in a liquid at room temperature with a time step of 2 fs, a skin of 0.1 nm with NL_STRIDE=10 is a safe choice.
  A larger skin allows a larger stride but makes every step slower, since more pairs are computed. You can check your choice by comparing a short run with one without neighbor list.

If GROUPB is empty, the coordination  number will be calculated based on the $\frac{N(N-1)}{2}$ pairs in GROUPA. This avoids computing
permuted indexes (e.g. pair (i,j) and (j,i)) twice and ensures that the calculation runs at twice the speed.

Notice that if there are common atoms between GROUPA and GROUPB the switching function should be
equal to one. These "self contacts" are discarded by plumed (since version 2.1),
so that they actually count as "zero".

## One atom against many atoms

To compute the coordination number of one atom (or a few) with many atoms, such as the oxygens of all the
water molecules around an ion, use COORDINATION with a neighbor list: between the updates of the list only
the atoms in the list are used, and only they are requested from the MD code.

```plumed
ow: GROUP ATOMS=2-30001:3
c: COORDINATION GROUPA=1 GROUPB=ow SWITCH={RATIONAL R_0=0.35 NN=6 MM=10 D_MAX=1.2} NLIST NL_CUTOFF=1.3 NL_STRIDE=10
```

If you need the coordination number of each atom of the small group as a vector, use [COORDINATIONNUMBER](COORDINATIONNUMBER.md)
with SPECIESA and SPECIESB.

## Examples

The following example instructs plumed to calculate the total coordination number of the atoms in group 1-10 with the atoms in group 20-100.  For atoms 1-10 coordination numbers are calculated that count the number of atoms from the second group that are within 0.3 nm of the central atom.  A neighbor list is used to make this calculation faster: it contains the pairs within 0.5 nm (a skin of 0.1 nm beyond D_MAX) and it is updated every 10 steps.

```plumed
COORDINATION GROUPA=1-10 GROUPB=20-100 SWITCH={RATIONAL R_0=0.3 D_MAX=0.4} NLIST NL_CUTOFF=0.5 NL_STRIDE=10
```

With larger groups the neighbor list is built faster with link cells:

```plumed
COORDINATION GROUPA=1-1000 GROUPB=1001-5000 SWITCH={RATIONAL R_0=0.3 D_MAX=0.4} NLISTCELLS NL_CUTOFF=0.5 NL_STRIDE=10
```

The following is a dummy example which should compute the value 0 because the self interaction
of atom 1 is skipped. Notice that in plumed 2.0 "self interactions" were not skipped, and the
same calculation should return 1.

```plumed
c: COORDINATION GROUPA=1 GROUPB=1 R_0=0.3
PRINT ARG=c STRIDE=10
````

Here's an example that shows what happens when COORDINATION is provided with
a single group:

```plumed
# define some huge group:
group: GROUP ATOMS=1-1000
# Here's coordination of a group against itself:
c1: COORDINATION GROUPA=group GROUPB=group R_0=0.3
# Here's coordination within a single group:
x: COORDINATION GROUPA=group R_0=0.3
# This is just multiplying times 2 the variable x:
c2: COMBINE ARG=x COEFFICIENTS=2 PERIODIC=NO

# the two variables c1 and c2 should be identical, but the calculation of c2 is twice faster
# since it runs on half of the pairs.
PRINT ARG=c1,c2 STRIDE=10
```

*/
//+ENDPLUMEDOC

class Coordination : public CoordinationBase {
  SwitchingFunction switchingFunction;

public:
  explicit Coordination(const ActionOptions&);
// active methods:
  static void registerKeywords( Keywords& keys );
  double pairing(double distance,double&dfunc,unsigned i,unsigned j)const override;
  double pairCutoff2() const override;
};

PLUMED_REGISTER_ACTION(Coordination,"COORDINATION")

void Coordination::registerKeywords( Keywords& keys ) {
  CoordinationBase::registerKeywords(keys);
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the switching functions that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.setValueDescription("scalar","the value of the coordination");
}

Coordination::Coordination(const ActionOptions&ao):
  Action(ao),
  CoordinationBase(ao) {

  std::string sw,errors;
  parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
    if( errors.length()!=0 ) {
      error("problem reading SWITCH keyword : " + errors );
    }
  } else {
    int nn=6;
    int mm=0;
    double d0=0.0;
    double r0=0.0;
    parse("R_0",r0);
    if(r0<=0.0) {
      error("R_0 should be explicitly specified and positive");
    }
    parse("D_0",d0);
    parse("NN",nn);
    parse("MM",mm);
    switchingFunction.set(nn,mm,r0,d0);
  }

  checkRead();

  log<<"  contacts are counted with cutoff "<<switchingFunction.description()<<"\n";
}

double Coordination::pairing(double distance,double&dfunc,unsigned i,unsigned j)const {
  (void) i; // avoid warnings
  (void) j; // avoid warnings
  return switchingFunction.calculateSqr(distance,dfunc);
}

double Coordination::pairCutoff2() const {
  // the switching functions are zero beyond dmax; the margin covers the rounding
  // of dmax*dmax and of the square root of the distance in calculateSqr
  return switchingFunction.get_dmax2()*(1.0+1e-12);
}

}

}
