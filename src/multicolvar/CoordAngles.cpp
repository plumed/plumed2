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

where $\theta_{ijk}$ is the angle between the vector connecting atom $i$ and $j$ and the vector connecting atom $j$ and atom $k$ and
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
private:
  /// Create the actions that calculate the angles between the bonds that connect catom to the atoms in group
  void createAngles( const std::string& lab, const std::string& catom, const std::vector<std::string>& group, const std::string& switch_input );
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

void CoordAngles::createAngles( const std::string& lab, const std::string& catom, const std::vector<std::string>& group, const std::string& switch_input ) {
  // Create the list of atoms.  An atom is not bonded to itself so the central atom is not used as one of its own neighbours
  std::string atlist;
  unsigned k=1;
  for(unsigned j=0; j<group.size(); ++j) {
    if( group[j]==catom ) {
      continue;
    }
    std::string num;
    Tools::convert( k, num );
    atlist += " ATOMS" + num + "=" + catom + "," + group[j];
    k++;
  }
  if( k==1 ) {
    error("there are no atoms in GROUP that are different from the central atom " + catom );
  }
  // Calculate the distances
  readInputLine( lab + "_dd: DISTANCE" + atlist );
  // Transform with the switching function
  readInputLine( lab + "_sw: LESS_THAN ARG=" + lab + "_dd SWITCH={" + switch_input +"}");
  // Now get the normalised vectors
  readInputLine( lab + "_comp: DISTANCE" + atlist + " COMPONENTS");
  readInputLine( lab + "_norm2: COMBINE ARG=" + lab + "_comp.x" + "," + lab + "_comp.y," + lab + "_comp.z POWERS=2,2,2 PERIODIC=NO");
  readInputLine( lab + "_norm: CUSTOM ARG=" + lab + "_norm2 FUNC=sqrt(x) PERIODIC=NO");
  readInputLine( lab + "_norm_x: CUSTOM ARG=" + lab + "_comp.x," + lab + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( lab + "_norm_y: CUSTOM ARG=" + lab + "_comp.y," + lab + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( lab + "_norm_z: CUSTOM ARG=" + lab + "_comp.z," + lab + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( lab + "_stack: VSTACK ARG=" + lab + "_norm_x" + "," + lab + "_norm_y," + lab + "_norm_z");
  readInputLine( lab + "_stackT: TRANSPOSE ARG=" + lab + "_stack");
  // Create the matrix of weights
  readInputLine( lab + "_swd: OUTER_PRODUCT ELEMENTS_ON_DIAGONAL_ARE_ZERO ARG=" + lab + "_sw," + lab + "_sw");
  // Avoid double counting
  readInputLine( lab + "_wmat: CUSTOM ARG=" + lab + "_swd FUNC=0.5*x PERIODIC=NO");
  // And the matrix of dot products and the angles
  readInputLine( lab + "_dpmat: MATRIX_PRODUCT ELEMENTS_ON_DIAGONAL_ARE_ZERO ARG=" + lab + "_stack," + lab + "_stackT");
  readInputLine( lab + "_angles: CUSTOM ARG=" + lab + "_dpmat FUNC=acos(x) PERIODIC=NO");
}

CoordAngles::CoordAngles(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Parse the central atoms
  std::vector<std::string> catoms;
  parseVector("CATOMS",catoms);
  Tools::interpretRanges(catoms);
  if( catoms.size()==0 ) {
    error("no central atoms were found in CATOMS");
  }
  // Parse the coordination sphere
  std::vector<std::string> group;
  parseVector("GROUP",group);
  Tools::interpretRanges(group);
  if( group.size()==0 ) {
    error("no atoms were found in GROUP");
  }
  std::string switch_input;
  parse("SWITCH",switch_input);
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
  const std::string lab = getShortcutLabel();
  if( catoms.size()==1 ) {
    createAngles( lab, catoms[0], group, switch_input );
    MultiColvarShortcuts::expandFunctions( lab, lab + "_angles", lab + "_wmat", keymap, this );
    if( do_mean ) {
      readInputLine( lab + "_denom: SUM ARG=" + lab + "_wmat PERIODIC=NO");
      readInputLine( lab + "_mean: CUSTOM ARG=" + lab + "_sum," + lab + "_denom FUNC=x/y PERIODIC=NO");
    }
    return;
  }
  // The angles are calculated separately for each of the central atoms, as the angles between bonds
  // that start from different atoms are not wanted.  Everything that we output is a sum over the angles
  // so the final quantities are the sums of the quantities for each central atom.
  std::vector<std::string> outputs, prefixes;
  for(unsigned i=0; i<catoms.size(); ++i) {
    std::string num;
    Tools::convert( i+1, num );
    prefixes.push_back( lab + "_c" + num );
    createAngles( prefixes[i], catoms[i], group, switch_input );
    std::vector<std::string> thisoutputs;
    MultiColvarShortcuts::expandFunctions( prefixes[i], prefixes[i] + "_angles", prefixes[i] + "_wmat", keymap, this, &thisoutputs );
    if( i==0 ) {
      outputs = thisoutputs;
    }
    if( do_mean ) {
      readInputLine( prefixes[i] + "_denom: SUM ARG=" + prefixes[i] + "_wmat PERIODIC=NO");
    }
  }
  for(const auto& suffix : outputs ) {
    std::string args = prefixes[0] + suffix;
    for(unsigned i=1; i<prefixes.size(); ++i) {
      args += "," + prefixes[i] + suffix;
    }
    readInputLine( lab + suffix + ": COMBINE ARG=" + args + " PERIODIC=NO");
  }
  if( do_mean ) {
    std::string args = prefixes[0] + "_denom";
    for(unsigned i=1; i<prefixes.size(); ++i) {
      args += "," + prefixes[i] + "_denom";
    }
    readInputLine( lab + "_denom: COMBINE ARG=" + args + " PERIODIC=NO");
    readInputLine( lab + "_mean: CUSTOM ARG=" + lab + "_sum," + lab + "_denom FUNC=x/y PERIODIC=NO");
  }
}

}
}
