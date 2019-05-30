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
#include "SymmetryFunctionBase.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace symfunc {

class AngularTetra : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit AngularTetra(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(AngularTetra,"TETRA_ANGULAR")

void AngularTetra::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::shortcutKeywords( keys ); 
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("compulsory","CUTOFF","-1","ignore distances that have a value larger than this cutoff");
  keys.remove("NN"); keys.remove("MM"); keys.remove("D_0"); keys.remove("R_0"); keys.remove("SWITCH");
}

AngularTetra::AngularTetra( const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  bool nopbc; parseFlag("NOPBC",nopbc); 
  std::string pbcstr=""; if( nopbc ) pbcstr = " NOPBC";
  // Read species input and create the matrix
  std::string sp_str, rcut; parse("SPECIES",sp_str); parse("CUTOFF",rcut);
  if( sp_str.length()>0 ) {
     readInputLine( getShortcutLabel() + "_mat: DISTANCE_MATRIX COMPONENTS GROUP=" + sp_str + " CUTOFF=" + rcut + pbcstr ); 
  } else {
     std::string specA, specB; parse("SPECIESA",specA); parse("SPECIESB",specB);
     if( specA.length()==0 ) error("missing input atoms");
     if( specB.length()==0 ) error("missing SPECIESB keyword");
     readInputLine( getShortcutLabel() + "_mat: DISTANCE_MATRIX COMPONENTS GROUPA=" + specA + " GROUPB=" + specB + " CUTOFF=" + rcut + pbcstr );
  }
  // Get the neighbors matrix
  readInputLine( getShortcutLabel() + "_neigh: NEIGHBORS ARG=" + getShortcutLabel() + "_mat.w NLOWEST=4");
  // Now construct the symmetry function (sum of cos(a) + 1/3)
  readInputLine( getShortcutLabel() + "_g8: GSYMFUNC_THREEBODY WEIGHT=" + getShortcutLabel() + "_neigh VECTORS1=" + getShortcutLabel() + "_mat.x " + 
                                      "VECTORS2=" + getShortcutLabel() + "_mat.y VECTORS3=" + getShortcutLabel() + "_mat.z FUNCTION1={TYPE=g8 LABEL=g8}");
  // Now evaluate the actual per atom CV
  readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_g8.g8 FUNC=(1-(3*x/8)) PERIODIC=NO"); 
  // And get the things to do with the quantities we have computed
  std::map<std::string,std::string> keymap; multicolvar::MultiColvarBase::readShortcutKeywords( keymap, this );
  multicolvar::MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

}
}
