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
#include "Torsion.h"
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR TORSION
/*
Calculate one or multiple torsional angles.

This command can be used to compute the torsion between four atoms as shown by the input below:

```plumed
t: TORSION ATOMS=1,2,3,4
PRINT ARG=t FILE=COLVAR
```

Alternatively you can use this action to calculate the angle between two vectors projected on the plane
orthogonal to an axis.  The example below uses this syntax and computes the cosine of the torsion that was calculated in the first example
input above.

```plumed
t: TORSION VECTORA=2,1 AXIS=2,3 VECTORB=3,4 COSINE
PRINT ARG=t FILE=COLVAR
```

If you combine this sytax with the functionality in [FIXEDATOM](FIXEDATOM.md) you can see how we can calculate the
torsional angle between two bond vectors around the z-axis as shown below:

```plumed
a0: FIXEDATOM AT=0,0,0
az: FIXEDATOM AT=0,0,1
t1: TORSION VECTORA=1,2 AXIS=a0,az VECTORB=5,6
PRINT ARG=t1 FILE=colvar STRIDE=20
```

If you are working with a protein you can specify the special named torsion angles $\phi$, $\psi$, $\omega$ and $\chi_1$
by using TORSION in combination with the [MOLINFO](MOLINFO.md) command.  This can be done by using the following
syntax.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=regtest/basic/rt32/helix.pdb
t1: TORSION ATOMS=@phi-3
t2: TORSION ATOMS=@psi-4
PRINT ARG=t1,t2 FILE=colvar STRIDE=10
```

Here, `@phi-3` tells plumed that you would like to calculate the $\phi$ angle in the third residue of the protein.
Similarly `@psi-4` tells plumed that you want to calculate the $\psi$ angle of the fourth residue of the protein.

If you would like to calculate multiple torsion angles at the same time you can use a command like the one shown below:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=regtest/basic/rt32/helix.pdb
t1: TORSION ATOMS1=@phi-3 ATOMS2=@phi-4 ATOMS3=@phi-5 ATOMS4=@phi-6 ATOMS5=@phi-7
PRINT ARG=t1 FILE=colvar STRIDE=10
```

This input tells PLUMED to calculate the $\phi$ angles in residues 3-7 of the protein.  The output, `t1`, is a 5 dimensional vector.

Notice that you can also use the VECTORA, VECTORB axis syntax when calculating multiple torsions as shown below:

```plumed
t: TORSION ...
  VECTORA1=2,1 AXIS1=2,3 VECTORB1=3,4
  VECTORA2=6,5 AXIS2=6,7 VECTORB2=7,8
  VECTORA3=10,9 AXIS3=10,11 VECTORB3=11,12
...
PRINT ARG=t FILE=colvar STRIDE=20
```

This input would output a three dimensional vector of torsion angles.

The last thing to note is that by default a procedure akin to that used in [WHOLEMOLECULES](WHOLEMOLECULES.md)
is used to ensure that the sets of atoms that are specified to each ATOMS keyword or set of VECTORA, AXIS and VECTORB keywords are not broken by the periodic
boundary conditions.  If you would like to turn this off for any reason you add the NOPBC in your input file as shown
below:

```plumed
t: TORSION ATOMS=1,2,3,4 NOPBC
PRINT ARG=t FILE=COLVAR
```


*/
//+ENDPLUMEDOC
typedef Torsion<double> TorsionD;
typedef ColvarShortcut<TorsionD> TorsionShortcut;
PLUMED_REGISTER_ACTION(TorsionShortcut,"TORSION")
PLUMED_REGISTER_ACTION(TorsionD,"TORSION_SCALAR")
typedef MultiColvarTemplate<TorsionD> TorsionMulti;
PLUMED_REGISTER_ACTION(TorsionMulti,"TORSION_VECTOR")

}
}

