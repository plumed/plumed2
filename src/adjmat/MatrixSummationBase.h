/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#ifndef __PLUMED_adjmat_MatrixSummationBase_h
#define __PLUMED_adjmat_MatrixSummationBase_h

#include "multicolvar/MultiColvarBase.h"
#include "AdjacencyMatrixVessel.h"

namespace PLMD {
namespace adjmat {

class MatrixSummationBase : public multicolvar::MultiColvarBase {
protected:
  /// The vessel that holds the adjacency matrix
  AdjacencyMatrixVessel* mymatrix;
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixSummationBase(const ActionOptions&);
  void calculate();
  void updateActiveAtoms( multicolvar::AtomValuePack& myatoms ) const ;
  bool isPeriodic(){ return false; }
  Vector getPositionOfAtomForLinkCells( const unsigned& iatom ) const ;
  bool isCurrentlyActive( const unsigned& bno, const unsigned& code );
};

}
}
#endif
