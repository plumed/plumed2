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
#include "tools/SwitchingFunction.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX DISTANCE_MATRIX
/*
Calculate a matrix of distances

\par Examples

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace adjmat {

class DistanceMatrix : public AdjacencyMatrixBase {
private:
  double cutoff;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DistanceMatrix(const ActionOptions&);
/// This does nothing
  double calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(DistanceMatrix,"DISTANCE_MATRIX")

void DistanceMatrix::registerKeywords( Keywords& keys ) {
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("compulsory","CUTOFF","-1","use a link cells algorithm with this cutoff to optimise the calculation - distances (but not derivatives) for atoms that are further apart than this cutoff and that in the same link cells will still be computed");
}

DistanceMatrix::DistanceMatrix( const ActionOptions& ao ):
  Action(ao),
  AdjacencyMatrixBase(ao) {
  // And set the link cell cutoff
  log.printf("  weight is distance between atoms \n");
  parse("CUTOFF",cutoff);
  if( cutoff<0 ) {
    setLinkCellCutoff( true, std::numeric_limits<double>::max() );
  } else {
    log.printf("  using link cells with cutoff %f to optimise calculation \n", cutoff);
    warning("some distances that are greater than the cutoff will still be evaluated. If you want to ignore them you will need to use a further CUSTOM action in your PLUMED input as detailed in the manual for version 2.11 or higher");
    setLinkCellCutoff( true, cutoff );
  }
}

double DistanceMatrix::calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const {
  Vector distance = pos2;
  double mod = distance.modulo();
  double invd = 1.0/mod;
  addAtomDerivatives( 0, (-invd)*distance, myvals );
  addAtomDerivatives( 1, (+invd)*distance, myvals );
  addBoxDerivatives( (-invd)*Tensor(distance,distance), myvals );
  return mod;
}

}
}

