/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/Action.h"
#include "core/ActionAnyorder.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Exception.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC INCLUDE
/*
Includes an external input file, similar to #include in C preprocessor.

The INCLUDE command is Useful when you want to split very large plumed.dat files.

In PLUMED 2.4 this action was now allowed to appear before the initial setup part of the file
(e.g. in the part that conained any [UNITS](UNITS.md) and [MOLINFO](MOLINFO.md) commands).
However, from PLUMED 2.5 onwards, INCLUDE commands can be used in any position of the file.

## Basic example

This input:

```plumed
c1: COM ATOMS=1-100
c2: COM ATOMS=101-202
d: DISTANCE ATOMS=c1,c2
PRINT ARG=d
```

can be replaced with this input:

```plumed
INCLUDE FILE=extras/pippo.dat
d: DISTANCE ATOMS=c1,c2
PRINT ARG=d
```

Notice that you can see the contents of pippo.dat by clicking on the name of the file

## Using INCLUDE to define groups

Using INCLUDE for the input in the previous section is rather pointless as the included file is rather short.
In a case like the one shown below the INCLUDE command is much more useful:

```plumed
INCLUDE FILE=extras/groups.dat
c: COORDINATION GROUPA=groupa GROUPB=groupb R_0=0.5
METAD ARG=c HEIGHT=0.2 PACE=100 SIGMA=0.2 BIASFACTOR=5
```

Although the `groups.dat` file here is short again it could be much larger if the groups
contained very lage numbers of atoms.

## INCLUDE and multiple replicas

Another case where INCLUDE can be useful is when running multi-replica simulations.
Here different replicas might have different input files, but perhaps a large part of the
input is shared. This part can be put in a common included file. For instance you could use
a `common.dat` to share part of the input across the input for the two replica simulation
which has the following `plumed.0.dat` input file:

```plumed
# this is plumed.0.dat
INCLUDE FILE=extras/common.dat
RESTRAINT ARG=t AT=1.0 KAPPA=10
```

and the following `plumed.1.dat` input file

```plumed
# this is plumed.1.dat
INCLUDE FILE=extras/common.dat
RESTRAINT ARG=t AT=1.2 KAPPA=10
```

When you run calculations in this way it is important to reember that PLUMED will always try to open
files for reading with a replica suffix first.  This is also even for files opened by INCLUDE!
We can thus implement the calculation above using the following common input for the two replicas:

```plumed
#SETTINGS NREPLICAS=2
t: TORSION ATOMS=1,2,3,4
INCLUDE FILE=extras/other.inc
```

We would then have an other.0.inc file that contains the defintion of the RESTRAINT command with AT=1.0
and an other.1.inc file that contains the defintion of the RESTRAINT command with AT=1.2.

In practice, however, we recommend using the special replica syntax that is discussed on [this page](parsing.md)
in place of the INCLUDE file approach that has been described above.  If you are using this syntax the input
for the calculation above is as follows:

```plumed
#SETTINGS NREPLICAS=2
t: TORSION ATOMS=1,2,3,4
RESTRAINT ARG=t AT=@replicas:1.0,1.2 KAPPA=10
```

The version of this input that uses INCLUDE files is an older syntax that was used in older versions of PLUMED
before the special replica syntax was made available.

*/
//+ENDPLUMEDOC

class Include :
  public ActionAnyorder {
public:
  static void registerKeywords( Keywords& keys );
  explicit Include(const ActionOptions&ao);
  void calculate() override {}
  void apply() override {}
};

PLUMED_REGISTER_ACTION(Include,"INCLUDE")

void Include::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  keys.add("compulsory","FILE","file to be included");
}

Include::Include(const ActionOptions&ao):
  Action(ao),
  ActionAnyorder(ao) {
  std::string f;
  parse("FILE",f);
  checkRead();
  plumed.readInputFile(f);
}

}
}

