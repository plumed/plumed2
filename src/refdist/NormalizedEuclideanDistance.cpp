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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC FUNCTION NORMALIZED_EUCLIDEAN_DISTANCE
/*
Calculate the normalised euclidean distance between two points in CV space

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace refdist {

class NormalizedEuclideanDistance : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit NormalizedEuclideanDistance(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(NormalizedEuclideanDistance,"NORMALIZED_EUCLIDEAN_DISTANCE")

void NormalizedEuclideanDistance::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG1","The poin that we are calculating the distance from");
  keys.add("compulsory","ARG2","The point that we are calculating the distance to");
  keys.add("compulsory","METRIC","The inverse covariance matrix that should be used when calculating the distance");
  keys.addFlag("SQUARED",false,"The squared distance should be calculated");
  keys.setValueDescription("the normalized euclidean distances between the input vectors");
  keys.needsAction("DISPLACEMENT"); keys.needsAction("CUSTOM"); keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("TRANSPOSE"); keys.needsAction("MATRIX_PRODUCT_DIAGONAL"); keys.needsAction("ONES");
}

NormalizedEuclideanDistance::NormalizedEuclideanDistance( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string arg1, arg2, metstr; parse("ARG1",arg1); parse("ARG2",arg2); parse("METRIC",metstr);
  // Vectors are in rows here
  readInputLine( getShortcutLabel() + "_diff: DISPLACEMENT ARG1=" + arg1 + " ARG2=" + arg2 );
  // Vectors are in columns here
  readInputLine( getShortcutLabel() + "_diffT: TRANSPOSE ARG=" + getShortcutLabel() + "_diff");
  // Get the action that computes the differences
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_diffT"); plumed_assert( av );
  // If this is a matrix we need create a matrix to multiply by
  if( av->copyOutput(0)->getRank()==2 ) {
    // Create some ones
    std::string nones; Tools::convert( av->copyOutput(0)->getShape()[1], nones );
    readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + nones);
    // Now do some multiplication to create a matrix that can be multiplied by our "inverse variance" vector
    if( av->copyOutput(0)->getShape()[0]==1 ) {
      readInputLine( getShortcutLabel() + "_" + metstr + "T: CUSTOM ARG=" + metstr + "," + getShortcutLabel() + "_ones FUNC=x*y PERIODIC=NO");
      readInputLine( getShortcutLabel() + "_" + metstr + ": TRANSPOSE ARG=" + getShortcutLabel() + "_" + metstr + "T");
    } else readInputLine( getShortcutLabel() + "_" + metstr + ": OUTER_PRODUCT ARG=" + metstr + "," + getShortcutLabel() + "_ones");
    metstr = getShortcutLabel() + "_" + metstr;
  }
  // Now do the multiplication
  readInputLine( getShortcutLabel() + "_sdiff: CUSTOM ARG=" + metstr + "," + getShortcutLabel() +"_diffT FUNC=x*y PERIODIC=NO");
  bool squared; parseFlag("SQUARED",squared); std::string olab = getShortcutLabel(); if( !squared ) olab += "_2";
  readInputLine( olab + ": MATRIX_PRODUCT_DIAGONAL ARG=" + getShortcutLabel() +"_diff," + getShortcutLabel() + "_sdiff");
  if( !squared ) readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_2 FUNC=sqrt(x) PERIODIC=NO");
}

}
}
