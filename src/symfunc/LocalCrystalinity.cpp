/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "CoordinationNumbers.h"
#include "multicolvar/MultiColvarShortcuts.h"

#include <complex>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR LOCAL_CRYSTALINITY
/*
Calculate the local crystalinity symmetry function

\par Examples


*/
//+ENDPLUMEDOC


class LocalCrystallinity : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit LocalCrystallinity(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(LocalCrystallinity,"LOCAL_CRYSTALINITY")

void LocalCrystallinity::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.add("numbered","GVECTOR","the coefficients of the linear combinations to compute for the CV");
  keys.setValueDescription("vector","the value of the local crystalinity for each of the input atoms");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
}

LocalCrystallinity::LocalCrystallinity( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // This builds an adjacency matrix
  std::string sp_str, specA, specB;
  parse("SPECIES",sp_str);
  parse("SPECIESA",specA);
  parse("SPECIESB",specB);
  CoordinationNumbers::expandMatrix( true, getShortcutLabel(), sp_str, specA, specB, this );
  // Input for denominator (coord)
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_mat.w," + getShortcutLabel() + "_ones");
  // Input for numerator
  std::string finput = "";
  for(unsigned i=1;; ++i) {
    std::vector<std::string> gvec;
    std::string istr;
    Tools::convert( i, istr );
    if( !parseNumberedVector("GVECTOR",i,gvec) ) {
      break;
    }
    if( gvec.size()!=3 ) {
      error("gvectors should have size 3");
    }
    // This is the dot product between the input gvector and the bond vector
    readInputLine( getShortcutLabel() + "_dot" + istr + ": COMBINE ARG=" + getShortcutLabel() + "_mat.x," + getShortcutLabel() + "_mat.y," + getShortcutLabel() + "_mat.z "
                   "PERIODIC=NO COEFFICIENTS=" + gvec[0] + "," + gvec[1] + "," + gvec[2] );
    // Now calculate the sine and cosine of the dot product
    readInputLine( getShortcutLabel() + "_cos" + istr + ": CUSTOM ARG=" + getShortcutLabel() +"_mat.w," + getShortcutLabel() + "_dot" + istr + " FUNC=x*cos(y) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_sin" + istr + ": CUSTOM ARG=" + getShortcutLabel() +"_mat.w," + getShortcutLabel() + "_dot" + istr + " FUNC=x*sin(y) PERIODIC=NO");
    // And sum up the sine and cosine over the coordination spheres
    readInputLine( getShortcutLabel() + "_cossum" + istr + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_cos" + istr + "," + getShortcutLabel() + "_ones");
    readInputLine( getShortcutLabel() + "_sinsum" + istr + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_sin" + istr + "," + getShortcutLabel() + "_ones");
    // And average the sine and cosine over the number of bonds
    readInputLine( getShortcutLabel() + "_cosmean" + istr + ": CUSTOM FUNC=x/y PERIODIC=NO ARG=" + getShortcutLabel() + "_cossum" + istr + "," + getShortcutLabel() + "_denom");
    readInputLine( getShortcutLabel() + "_sinmean" + istr + ": CUSTOM FUNC=x/y PERIODIC=NO ARG=" + getShortcutLabel() + "_sinsum" + istr + "," + getShortcutLabel() + "_denom");
    // And work out the square modulus of this complex number
    readInputLine( getShortcutLabel() + "_" + istr + ": CUSTOM FUNC=x*x+y*y PERIODIC=NO ARG=" + getShortcutLabel() + "_cosmean" + istr + "," + getShortcutLabel() + "_sinmean" + istr);
    // These are all the kvectors that we are adding together in the final combine for the final CV
    if( i>1 ) {
      finput += ",";
    }
    finput += getShortcutLabel() + "_" + istr;
  }
  // This computes the final CV
  readInputLine( getShortcutLabel() + ": COMBINE NORMALIZE PERIODIC=NO ARG=" + finput );
  // Now calculate the total length of the vector
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", this );
}

}
}

