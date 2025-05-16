/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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

//+PLUMEDOC MCOLVAR VORONOI
/*
Do a voronoi analysis

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class Voronoi : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit Voronoi(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Voronoi,"VORONOI")

void Voronoi::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","matrix","the distance/adjacency matrix that should be used to perform the voronoi analysis");
  keys.setValueDescription("matrix","a matrix in which element ij is equal to one if the ij component of the input matrix is lower than all the ik elements of the matrix where k is not j and zero otherwise");
  keys.needsAction("NEIGHBORS");
}

Voronoi::Voronoi(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  readInputLine( getShortcutLabel() + ": NEIGHBORS NLOWEST=1 " + convertInputLineToString() );
}

}
}
