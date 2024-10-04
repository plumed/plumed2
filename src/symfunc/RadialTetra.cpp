/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/ActionWithValue.h"
#include "CoordinationNumbers.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC MCOLVAR TETRA_RADIAL
/*
Calculate the radial tetra CV

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace symfunc {

class RadialTetra : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit RadialTetra(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(RadialTetra,"TETRA_RADIAL")

void RadialTetra::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("compulsory","CUTOFF","-1","ignore distances that have a value larger than this cutoff");
  keys.remove("NN");
  keys.remove("MM");
  keys.remove("D_0");
  keys.remove("R_0");
  keys.remove("SWITCH");
  keys.needsAction("DISTANCE_MATRIX");
  keys.needsAction("NEIGHBORS");
  keys.needsAction("CUSTOM");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
}

RadialTetra::RadialTetra( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read species input and create the matrix
  bool nopbc;
  parseFlag("NOPBC",nopbc);
  std::string pbcstr="";
  if( nopbc ) {
    pbcstr = " NOPBC";
  }
  std::string sp_str, rcut;
  parse("SPECIES",sp_str);
  parse("CUTOFF",rcut);
  if( sp_str.length()>0 ) {
    readInputLine( getShortcutLabel() + "_mat: DISTANCE_MATRIX GROUP=" + sp_str + " CUTOFF=" + rcut + pbcstr );
  } else {
    std::string specA, specB;
    parse("SPECIESA",specA);
    parse("SPECIESB",specB);
    if( specA.length()==0 ) {
      error("missing input atoms");
    }
    if( specB.length()==0 ) {
      error("missing SPECIESB keyword");
    }
    readInputLine( getShortcutLabel() + "_mat: DISTANCE_MATRIX GROUPA=" + specA + " GROUPB=" + specB + " CUTOFF=" + rcut + pbcstr);
  }
  // Get the neighbors matrix
  readInputLine( getShortcutLabel() + "_neigh: NEIGHBORS ARG=" + getShortcutLabel() + "_mat NLOWEST=4");
  // Now get distance matrix that just contains four nearest distances
  readInputLine( getShortcutLabel() + "_near4: CUSTOM ARG=" + getShortcutLabel() + "_mat," + getShortcutLabel() + "_neigh FUNC=x*y PERIODIC=NO");
  //Now compute sum of four nearest distances
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_sum4: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_near4," + getShortcutLabel() + "_ones");
  // Now compute squares of four nearest distance
  readInputLine( getShortcutLabel() + "_near4_2: CUSTOM ARG=" + getShortcutLabel() + "_near4 FUNC=x*x PERIODIC=NO");
  // Now compute sum of the squares of the four nearest distances
  readInputLine( getShortcutLabel() + "_sum4_2: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_near4_2," + getShortcutLabel() + "_ones");
  // Evaluate the average distance to the four nearest neighbors
  readInputLine( getShortcutLabel() + "_meanr: CUSTOM ARG=" + getShortcutLabel() + "_sum4 FUNC=0.25*x PERIODIC=NO");
  // Now evaluate the actual per atom CV
  readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_sum4," + getShortcutLabel() + "_sum4_2," + getShortcutLabel() + "_meanr " +
                 "FUNC=(1-(y-x*z)/(12*z*z)) PERIODIC=NO");
  // And get the things to do with the quantities we have computed
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

}
}
