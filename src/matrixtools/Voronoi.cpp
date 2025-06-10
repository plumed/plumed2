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

This shortcut allows you to perform a [Voronoi analysis](https://en.wikipedia.org/wiki/Voronoi_diagram).
The input to this action is a rectangular matrix that describes the distances between two sets of points.
The points for which this action is receiving distances between could be atom positions or they could represent the
dissimilarities between various trajectory frames that have been stored using [COLLECT_FRAMES](COLLECT_FRAMES.md).
In the example below the matrix of input distances are distances between the positions of atoms:

```plumed
d: DISTANCE_MATRIX GROUPA=1-10 GROUPB=11-100
v: VORONOI ARG=d
ones: ONES SIZE=90
c: MATRIX_VECTOR_PRODUCT ARG=v,ones
PRINT ARG=c FILE=colvar
```

The VORONOI action outputs a rectangular matrix, $V$, with the same shape as the input matrix, $D$. $V_{ij}$
is 1 if $D_{ij}$ is the shortest distance in row $i$.


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
