/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "DihedralCorrelation.h"
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "tools/Torsion.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DIHEDRAL_CORRELATION
/*
Measure the correlation between a pair of dihedral angles

This CV measures the correlation between two dihedral angles, $\phi$ and $\psi$, as follows:

$$
s = \frac{1}{2} \left[ 1 + \cos( \phi - \psi ) \right]
$$

An example input that computes and calculates this quantity with $\phi$ being the dihedral
angle involving atoms 1, 2, 3 and 4 and $\psi$ being the angle involving atoms 7, 8, 9 and 10
is shown below:

```plumed
d: DIHEDRAL_CORRELATION ATOMS=1,2,3,4,5,6,7,8
PRINT ARG=d FILE=colvar
```

If you want to calculate the dihedral correlations between multiple pairs of dihedral angles using
this action you would use an input like this one shown below:

```plumed
d: DIHEDRAL_CORRELATION ...
  ATOMS1=1,2,3,4,5,6,7,8
  ATOMS2=9,10,11,12,13,14,15,16
...
PRINT ARG=d FILE=colvar
```

This input calculates and outputs a two dimensional vector that contains two of these dihedral correlation
values. Commands similar to these are used within the [DIHCOR](DIHCOR.md) shortcut.

The last thing to note is that by default a procedure akin to that used in [WHOLEMOLECULES](WHOLEMOLECULES.md)
is used to ensure that the sets of atoms that are specified to each ATOMS keyword are not broken by the periodic
boundary conditions.  If you would like to turn this off for any reason you add the NOPBC in your input file as shown
below:

```plumed
d: DIHEDRAL_CORRELATION ATOMS=1,2,3,4,5,6,7,8 NOPBC
PRINT ARG=d FILE=colvar
```

*/
//+ENDPLUMEDOC


typedef DihedralCorrelation<double> DihedralCorrelationD;
typedef ColvarShortcut<DihedralCorrelationD> DihedralCorrelationShortcut;
PLUMED_REGISTER_ACTION(DihedralCorrelationShortcut,"DIHEDRAL_CORRELATION")
PLUMED_REGISTER_ACTION(DihedralCorrelationD,"DIHEDRAL_CORRELATION_SCALAR")
typedef MultiColvarTemplate<DihedralCorrelationD> DihedralCorrelationMulti;
PLUMED_REGISTER_ACTION(DihedralCorrelationMulti,"DIHEDRAL_CORRELATION_VECTOR")

}
}
