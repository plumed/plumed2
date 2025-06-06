/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "MultiColvarShortcuts.h"

//+PLUMEDOC MCOLVAR COORD_ANGLES
/*
Calculate functions of the distribution of angles between bonds in the first coordination spheres of a set of atoms

This action ncan be used to calculate functions, $g$, such as:

$$
f(x) = \sum_{ijk} s(r_{ij})s(r_{jk}) g(\theta_{ijk})
$$

where $\theta_{ijk}$ is the angle between the vector connecting atom $i$ and and $j$ and the vector connecting atom $j$ and atom $k$ and
where $s(r)$ is a switching function.  The switching functions in the expression above ensure that we can calculate all the angles in the first coordination
spheres of an atom using an input like the one shown below:

```plumed
c1: COORD_ANGLES ...
  CATOMS=1 GROUP=2-100 SWITCH={RATIONAL R_0=1.0} SUM
...
PRINT ARG=c1.sum FILE=colvar
```

The input above will output the sum of all the angles in the first coordination sphere.

__As you can see if you expand the shortcut in the input above, the calculation of the above sum is obtained by calculating a [Hadamard product](https://en.wikipedia.org/wiki/Hadamard_product_(matrices))
of an [OUTER_PRODUCT](OUTER_PRODUCT.md) matrix and a [MATRIX_PRODUCT](MATRIX_PRODUCT.md).  We have thus deprecated this action as we believe that there is value in learning to use the more complicated syntax
directly.__

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class CoordAngles : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  static void pruneShortcuts(Keywords& keys);
  explicit CoordAngles(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(CoordAngles,"COORD_ANGLES")

void CoordAngles::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","CATOMS","all the angles between the bonds that radiate out from these central atom are computed");
  keys.add("atoms","GROUP","a list of angls between pairs of bonds connecting one of the atoms specified using the CATOM command and two of the atoms specified here are computed");
  keys.add("compulsory","SWITCH","the switching function specifies that only those bonds that have a length that is less than a certain threshold are considered");
  MultiColvarShortcuts::shortcutKeywords( keys );
  pruneShortcuts( keys );
  keys.needsAction("DISTANCE");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
  keys.needsAction("VSTACK");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("MATRIX_PRODUCT");
  keys.setDeprecated("MATRIX_PRODUCT");
}

void CoordAngles::pruneShortcuts(Keywords& keys) {
  keys.remove("ALT_MIN");
  keys.remove("MIN");
  keys.remove("MAX");
  keys.remove("HIGHEST");
  keys.remove("LOWEST");
}

CoordAngles::CoordAngles(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Parse the central atoms
  std::vector<std::string> catoms;
  parseVector("CATOMS",catoms);
  Tools::interpretRanges(catoms);
  // Parse the coordination sphere
  std::vector<std::string> group;
  parseVector("GROUP",group);
  Tools::interpretRanges(group);
  // Create the list of atoms
  std::string atlist;
  unsigned k=1;
  for(unsigned i=0; i<catoms.size(); ++i) {
    for(unsigned j=0; j<group.size(); ++j) {
      std::string num;
      Tools::convert( k, num );
      atlist += " ATOMS" + num + "=" + catoms[i] + "," + group[j];
      k++;
    }
  }
  // Calculate the distances
  readInputLine( getShortcutLabel() + "_dd: DISTANCE" + atlist );
  // Transform with the switching function
  std::string switch_input;
  parse("SWITCH",switch_input);
  readInputLine( getShortcutLabel() + "_sw: LESS_THAN ARG=" + getShortcutLabel() + "_dd SWITCH={" + switch_input +"}");
  // Now get the normalised vectors
  readInputLine( getShortcutLabel() + "_comp: DISTANCE" + atlist + " COMPONENTS");
  readInputLine( getShortcutLabel() + "_norm2: COMBINE ARG=" + getShortcutLabel() + "_comp.x" + "," + getShortcutLabel() + "_comp.y," + getShortcutLabel() + "_comp.z POWERS=2,2,2 PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm: CUSTOM ARG=" + getShortcutLabel() + "_norm2 FUNC=sqrt(x) PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm_x: CUSTOM ARG=" + getShortcutLabel() + "_comp.x," + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm_y: CUSTOM ARG=" + getShortcutLabel() + "_comp.y," + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm_z: CUSTOM ARG=" + getShortcutLabel() + "_comp.z," + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_stack: VSTACK ARG=" + getShortcutLabel() + "_norm_x" + "," + getShortcutLabel() + "_norm_y," + getShortcutLabel() + "_norm_z");
  readInputLine( getShortcutLabel() + "_stackT: TRANSPOSE ARG=" + getShortcutLabel() + "_stack");
  // Create the matrix of weights
  readInputLine( getShortcutLabel() + "_swd: OUTER_PRODUCT ELEMENTS_ON_DIAGONAL_ARE_ZERO ARG=" + getShortcutLabel() + "_sw," + getShortcutLabel() + "_sw");
  // Avoid double counting
  readInputLine( getShortcutLabel() + "_wmat: CUSTOM ARG=" + getShortcutLabel() + "_swd FUNC=0.5*x PERIODIC=NO");
  // And the matrix of dot products and the angles
  readInputLine( getShortcutLabel() + "_dpmat: MATRIX_PRODUCT ELEMENTS_ON_DIAGONAL_ARE_ZERO ARG=" + getShortcutLabel() + "_stack," + getShortcutLabel() + "_stackT");
  readInputLine( getShortcutLabel() + "_angles: CUSTOM ARG=" + getShortcutLabel() + "_dpmat FUNC=acos(x) PERIODIC=NO");
  // Read the input
  Keywords keys;
  MultiColvarShortcuts::shortcutKeywords( keys );
  pruneShortcuts( keys );
  bool do_mean;
  parseFlag("MEAN",do_mean);
  std::map<std::string,std::string> keymap;
  readShortcutKeywords( keys, keymap );
  if( do_mean ) {
    keymap.insert(std::pair<std::string,std::string>("SUM",""));
  }
  MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_angles", getShortcutLabel() + "_wmat", keymap, this );
  if( do_mean ) {
    readInputLine( getShortcutLabel() + "_denom: SUM ARG=" + getShortcutLabel() + "_wmat PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_mean: CUSTOM ARG=" + getShortcutLabel() + "_sum," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  }
}

}
}
