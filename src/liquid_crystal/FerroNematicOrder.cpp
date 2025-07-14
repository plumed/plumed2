/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 of Alexander Humeniuk.

   This file is part of the liquid_crystal plumed module.

   The liquid_crystal plumed module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The liquid_crystal plumed module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "multicolvar/MultiColvarShortcuts.h"

using namespace PLMD::multicolvar;

namespace PLMD {
namespace liquid_crystal {

//+PLUMEDOC COLVAR FERRONEMATIC_ORDER
/*
Calculate the ferronematic order parameter.

The ferronematic order parameter P depends on the relative orientations of the molecular
axes. If the axes all point into the same direction, giving rise to a net polarization if
the molecules have permanent dipole moments, P is close to 1. If the molecular axes are
oriented isotropically or are aligned antiparallel so that there is no net polarization,
P is close 0.

The nematic and ferronematic order parameters can be used to distinguish the isotropic,
nematic and ferronematic phases of liquid crystals.

$P$ is length of the average of the molecular axes ($\hat{u}_i$ for $i=1,\ldots,N$)
$$
P = \vert \frac{1}{N} \sum_{i=1}^N \hat{u}_i \vert
$$
Since the molecular axes are unit vectors, P ranges from a minimum of 0 to a maximum of 1.

By adding a bias to the ferronematic order parameter, one can drive a liquid crystal from the
isotropic to the ferronematic phase.

The axis of a rod-like molecule is defined as the distance vector between two atoms,
it points from the tail atom to the head atom.

```plumed
# Assume there are three molecules with 20 atoms each.
# In the first molecule the molecular axis vector points from atom 1 to atom 20,
# in the second molecule it points from atom 21 to atom 40
# and in the third from atom 41 to atom 60.
# The ferronematic order parameter for the three molecules is computed as
P: FERRONEMATIC_ORDER MOLECULE_STARTS=1,21,41 MOLECULE_ENDS=20,40,60
PRINT FILE=colvar ARG=P

# Add a bias to the ferronematic order parameter P.
BIASVALUE ARG=P
```

*/
//+ENDPLUMEDOC

class FerroNematicOrder : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit FerroNematicOrder(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(FerroNematicOrder,"FERRONEMATIC_ORDER")

void FerroNematicOrder::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","MOLECULE_STARTS","The atoms where the molecular axis starts.");
  keys.add("atoms","MOLECULE_ENDS","The atoms where the molecular axis ends.");
  keys.setValueDescription("scalar","the modulus of the average vector");
  keys.needsAction("DISTANCE");
  keys.needsAction("CUSTOM");
  keys.needsAction("MEAN");
}

FerroNematicOrder:: FerroNematicOrder(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Fetch indices of atoms that define the tails and the heads of the molecular axes.
  std::vector<std::string> starts, ends;
  MultiColvarShortcuts::parseAtomList("MOLECULE_STARTS",starts,this);
  MultiColvarShortcuts::parseAtomList("MOLECULE_ENDS",ends,this);

  if( starts.size()!=ends.size() )
    error(
      "Mismatched numbers of atoms specified to MOLECULE_STARTS and MOLECULE_ENDS keywords. "
      "The molecular axes are specified by pairs of atoms."
    );

  std::string dlist = "";
  for(unsigned i=0; i<starts.size(); ++i) {
    std::string num;
    Tools::convert( i+1, num );
    dlist += " ATOMS" + num + "=" + starts[i] + "," + ends[i];
  }

  std::string L = getShortcutLabel();
  // Calculate the lengths of the distance vectors
  //   d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ...
  readInputLine( L + "_dvals: DISTANCE" + dlist );
  // Calculate the molecular axes of the molecules
  //   dc: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=3,4 ...
  readInputLine( L + "_dvecs: DISTANCE COMPONENTS " + dlist );
  // Convert the molecular axes into unit vectors
  //   dux: CUSTOM ARG=dc.x,d FUNC=x/y PERIODIC=NO
  //   duy: CUSTOM ARG=dc.y,d FUNC=x/y PERIODIC=NO
  //   duz: CUSTOM ARG=dc.z,d FUNC=x/y PERIODIC=NO
  readInputLine( L + "_dux: CUSTOM ARG=" + L + "_dvecs.x," + L + "_dvals FUNC=x/y PERIODIC=NO");
  readInputLine( L + "_duy: CUSTOM ARG=" + L + "_dvecs.y," + L + "_dvals FUNC=x/y PERIODIC=NO");
  readInputLine( L + "_duz: CUSTOM ARG=" + L + "_dvecs.z," + L + "_dvals FUNC=x/y PERIODIC=NO");
  // Now calculate the average of the molecular axes
  //   mux: MEAN ARG=dux PERIODIC=NO
  //   muy: MEAN ARG=duz PERIODIC=NO
  //   muz: MEAN ARG=dyz PERIODIC=NO
  readInputLine( L + "_mux: MEAN ARG=" + L + "_dux PERIODIC=NO");
  readInputLine( L + "_muy: MEAN ARG=" + L + "_duy PERIODIC=NO");
  readInputLine( L + "_muz: MEAN ARG=" + L + "_duz PERIODIC=NO");
  // Compute the ferronematic order parameter
  //   p: CUSTOM ARG=mux,muy,muz FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
  readInputLine( L + ": CUSTOM ARG=" + L + "_mux," + L + "_muy," + L + "_muz FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO");
}

}
}
