/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "AdjacencyMatrixBase.h"
#include "core/ActionRegister.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX DISTANCE_MATRIX
/*
Calculate a matrix of distances

\par Examples

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace adjmat {

class DistanceMatrix {
public:
  double cutoff;
  static void registerKeywords( Keywords& keys );
  void parseInput( AdjacencyMatrixBase<DistanceMatrix>* action );
  DistanceMatrix & operator=( const DistanceMatrix& m ) {
    cutoff = m.cutoff;
    return *this;
  }
  static void calculateWeight( const DistanceMatrix& data, const AdjacencyMatrixInput& input, MatrixOutput& output );
};

typedef AdjacencyMatrixBase<DistanceMatrix> dmap;
PLUMED_REGISTER_ACTION(dmap,"DISTANCE_MATRIX")

void DistanceMatrix::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","CUTOFF","-1","ignore distances that have a value larger than this cutoff");
}

void DistanceMatrix::parseInput( AdjacencyMatrixBase<DistanceMatrix>* action ) {
  // And set the link cell cutoff
  action->log.printf("  weight is distance between atoms \n");
  action->parse("CUTOFF",cutoff);
  if( cutoff<0 ) {
    action->setLinkCellCutoff( true, std::numeric_limits<double>::max() );
  } else {
    action->log.printf("  ignoring distances that are larger than %f \n", cutoff);
    action->setLinkCellCutoff( true, cutoff );
  }
}

void DistanceMatrix::calculateWeight( const DistanceMatrix& data, const AdjacencyMatrixInput& input, MatrixOutput& output ) {
  output.val[0] = input.pos.modulo();
  if( data.cutoff<0 || output.val[0]<data.cutoff ) {
    double invd = 1.0/output.val[0];
    Vector v = (-invd)*input.pos;
    output.deriv[0] = v[0];
    output.deriv[1] = v[1];
    output.deriv[2] = v[2];
    output.deriv[3] = -v[0];
    output.deriv[4] = -v[1];
    output.deriv[5] = -v[2];
    Tensor t = (-invd)*Tensor(input.pos,input.pos);
    output.deriv[6] = t[0][0];
    output.deriv[7] = t[0][1];
    output.deriv[8] = t[0][2];
    output.deriv[9] = t[1][0];
    output.deriv[10] = t[1][1];
    output.deriv[11] = t[1][2];
    output.deriv[12] = t[2][0];
    output.deriv[13] = t[2][1];
    output.deriv[14] = t[2][2];
  }
}

}
}

