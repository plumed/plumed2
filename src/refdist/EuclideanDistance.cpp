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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"

//+PLUMEDOC MCOLVAR EUCLIDEAN_DISTANCE
/*
Calculate the euclidean distance between two vectors of arguments

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace refdist {

class EuclideanDistance : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit EuclideanDistance(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(EuclideanDistance,"EUCLIDEAN_DISTANCE")

void EuclideanDistance::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG1","The poin that we are calculating the distance from");
  keys.add("compulsory","ARG2","The point that we are calculating the distance to");
  keys.addFlag("SQUARED",false,"The squared distance should be calculated");
  keys.setValueDescription("scalar/vector","the euclidean distances between the input vectors");
  keys.needsAction("DISPLACEMENT");
  keys.needsAction("CUSTOM");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("MATRIX_PRODUCT_DIAGONAL");
}

EuclideanDistance::EuclideanDistance( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string arg1, arg2;
  parse("ARG1",arg1);
  parse("ARG2",arg2);
  // Vectors are in rows here
  readInputLine( getShortcutLabel() + "_diff: DISPLACEMENT ARG1=" + arg1 + " ARG2=" + arg2 );
  // Get the action that computes the differences
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_diff");
  plumed_assert( av );
  // Check if squared
  bool squared;
  parseFlag("SQUARED",squared);
  std::string olab = getShortcutLabel();
  if( !squared ) {
    olab += "_2";
  }
  // Deal with an annoying corner case when displacement has a single argument
  if( av->copyOutput(0)->getRank()==0 ) {
    readInputLine( olab + ": CUSTOM ARG=" + getShortcutLabel() + "_diff FUNC=x*x PERIODIC=NO");
  } else {
    // Notice that the vectors are in the columns here
    readInputLine( getShortcutLabel() + "_diffT: TRANSPOSE ARG=" + getShortcutLabel() + "_diff");
    readInputLine( olab + ": MATRIX_PRODUCT_DIAGONAL ARG=" + getShortcutLabel() + "_diff," + getShortcutLabel() + "_diffT");
  }
  if( !squared ) {
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_2 FUNC=sqrt(x) PERIODIC=NO");
  }
}

}
}
