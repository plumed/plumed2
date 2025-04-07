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
#include "PathMSDBase.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR PATHMSD
/*
This Colvar calculates path collective variables.

Path collective variables were introduced in the paper cited in the bibliography.
This CV defines a a curvilinear path using a set of reference configurations. Variables
that measure the instantaneous position, sss, of the system on this path and the distance from
this path, zzz, are then computed.  The following input illustrates the syntax that is used
with this variable

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt39/all.pdb
p1: PATHMSD REFERENCE=regtest/basic/rt39/all.pdb LAMBDA=500.0 NEIGH_STRIDE=4 NEIGH_SIZE=8
PRINT ARG=p1.sss,p1.zzz STRIDE=1 FILE=colvar
```

The `NEIGH_STRIDE=4` and `NEIGH_SIZE=8` keywords here control the neighbor list parameter (optional but
recommended for performance) and states that the neighbor list will be calculated every 4
steps and consider only the closest 8 member to the actual md snapshots.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding molecules with a procedure
that is equivalent to that done in [WHOLEMOLECULES](WHOLEMOLECULES.md). Notice that
rebuilding is local to this action. This is different from [WHOLEMOLECULES](WHOLEMOLECULES.md)
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

The implementation of this collective variable and of [PROPERTYMAP](PROPERTYMAP.md)
is shared, as well as most input options.


*/
//+ENDPLUMEDOC

class PathMSD : public PathMSDBase {
public:
  explicit PathMSD(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(PathMSD,"PATHMSD")

void PathMSD::registerKeywords(Keywords& keys) {
  PathMSDBase::registerKeywords(keys);
  keys.addOutputComponent("sss","default","scalar","the position on the path");
  keys.addOutputComponent("zzz","default","scalar","the distance from the path");
  keys.addDOI("10.1063/1.2432340");
}

PathMSD::PathMSD(const ActionOptions&ao):
  Action(ao),PathMSDBase(ao) {
  checkRead();

  log<<"  Bibliography "
     <<plumed.cite("Branduardi, Gervasio, Parrinello J. Chem. Phys. 126, 054103 (2007)")
     <<"\n";
  // no need to read anything
  addComponentWithDerivatives("sss");
  componentIsNotPeriodic("sss");
  addComponentWithDerivatives("zzz");
  componentIsNotPeriodic("zzz");
  requestAtoms(pdbv[0].getAtomNumbers());

  double i=1.;
  for(unsigned it=0 ; it<nframes ; ++it) {
    std::vector<double> v;
    v.push_back(i);
    indexvec.push_back(v);
    i+=1.;
  }
}

}

}
