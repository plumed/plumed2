/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "MultiColvarBase.h"
#include <string>
#include <cmath>

namespace PLMD {
namespace multicolvar {

class AlphaBeta : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit AlphaBeta(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(AlphaBeta,"ALPHABETA")

void AlphaBeta::registerKeywords(Keywords& keys) {
  MultiColvarBase::registerKeywords( keys ); 
  keys.add("compulsory","REFERENCE","the reference values for each of the torsional angles.  If you use a single REFERENCE value the "
           "same reference value is used for all torsions");
}

AlphaBeta::AlphaBeta(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Read in the reference value
  std::string refstr; parse("REFERENCE",refstr);
  // Calculate angles
  readInputLine( getShortcutLabel() + "_torsions: TORSIONS " + convertInputLineToString() ); 
  // Caculate difference from reference using combine
  readInputLine( getShortcutLabel() + "_comb: COMBINE PARAMETERS=" + refstr + " ARG1=" + getShortcutLabel() + "_torsions PERIODIC=NO" );
  // Now matheval for cosine bit
  readInputLine( getShortcutLabel() + "_cos: MATHEVAL ARG1=" + getShortcutLabel() + "_comb FUNC=0.5+0.5*cos(x) PERIODIC=NO");
  // And combine to get final value
  readInputLine( getShortcutLabel() + ": SUM ARG=" + getShortcutLabel() + "_cos PERIODIC=NO");
}

}
}
