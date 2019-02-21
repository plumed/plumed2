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

namespace PLMD {
namespace function {

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
}

NormalizedEuclideanDistance::NormalizedEuclideanDistance( const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{ 
  std::string arg1, arg2, metstr; parse("ARG1",arg1); parse("ARG2",arg2); parse("METRIC",metstr);
  readInputLine( getShortcutLabel() + "_diff: DIFFERENCE ARG1=" + arg1 + " ARG2=" + arg2 );
  readInputLine( getShortcutLabel() + "_vdot: MATHEVAL ARG1=" + metstr + " ARG2=" + getShortcutLabel() +"_diff FUNC=x*y*y PERIODIC=NO");
  bool squared; parseFlag("SQUARED",squared);
  if( !squared ) {
      readInputLine( getShortcutLabel() + "_2: COMBINE ARG=" + getShortcutLabel() + "_vdot PERIODIC=NO");
      readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_2 FUNC=sqrt(x) PERIODIC=NO");
  } else readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_vdot PERIODIC=NO");
}

}
}
