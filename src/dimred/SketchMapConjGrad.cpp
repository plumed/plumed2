/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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

//+PLUMEDOC DIMRED SKETCHMAP_CONJGRAD
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class SketchMapConjGrad : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit SketchMapConjGrad( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(SketchMapConjGrad,"SKETCHMAP_CONJGRAD")

void SketchMapConjGrad::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("numbered","ARG","the starting coordinate for each projection");
  keys.add("compulsory","WEIGHTS","a vector containing the weights of the points");
  keys.add("compulsory","DISSIMILARITIES","the matrix of dissimilarities that are to be reproduced");
  keys.add("compulsory","HIGH_DIM_FUNCTION","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","the parameters of the switching function in the low dimensional space");
}

SketchMapConjGrad::SketchMapConjGrad( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  // Transform the dissimilarities using the switching function
  std::string dissmat, hdfunc; parse("DISSIMILARITIES",dissmat); parse("HIGH_DIM_FUNCTION",hdfunc);
  readInputLine( getShortcutLabel() + "_hdmat: MORE_THAN ARG1=" + dissmat + " SQUARED SWITCH={" + hdfunc + "}");
  // Now for the weights - read the vector of weights first
  std::string wvec; parse("WEIGHTS",wvec);
  // Now calculate the sum of thse weights
  readInputLine( wvec + "_sum: COMBINE ARG=" + wvec + " PERIODIC=NO");
  // And normalise the vector of weights using this sum
  readInputLine( wvec + "_normed: CUSTOM ARG1=" + wvec + "_sum ARG2=" + wvec + " FUNC=y/x PERIODIC=NO");
  // And now create the matrix of weights
  readInputLine( wvec + "_mat: DOTPRODUCT_MATRIX GROUP1=" + wvec + "_normed");
  // Run the arrange points object
  std::string ldfunc; parse("LOW_DIM_FUNCTION",ldfunc);
  readInputLine( getShortcutLabel() + ": ARRANGE_POINTS " + convertInputLineToString() + " TARGET1=" + getShortcutLabel() + "_hdmat FUNC1={" + ldfunc + "} WEIGHTS1=" + wvec + "_mat");  
}

}
}
