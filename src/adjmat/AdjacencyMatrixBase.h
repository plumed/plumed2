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
#ifndef __PLUMED_adjmat_AdjacencyMatrixBase_h
#define __PLUMED_adjmat_AdjacencyMatrixBase_h

#include <vector>
#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "tools/LinkCells.h"
#include "MatrixElementPack.h"

namespace PLMD {
namespace adjmat {

class AdjacencyMatrixBase : 
  public ActionAtomistic,
  public ActionWithValue 
{
friend class MatrixElementPack;
private:
  bool nopbc, components;
  LinkCells linkcells;
  std::vector<unsigned> ablocks;
  void updateWeightDerivativeIndices( const unsigned& sno, const std::vector<unsigned>& indices, MultiValue& myvals ) const ;
protected:
  void setLinkCellCutoff( const double& lcut, double tcut=-1.0 );
public:
  static void registerKeywords( Keywords& keys );
  explicit AdjacencyMatrixBase(const ActionOptions&);
  unsigned getNumberOfDerivatives() const ;
  void buildCurrentTaskList( std::vector<unsigned>& tflags );
  void calculate();
  void performTask( const unsigned& task_index, MultiValue& myvals ) const ;
  virtual double calculateWeight( MatrixElementPack& myvals ) const = 0;
  void apply();
};

inline
unsigned AdjacencyMatrixBase::getNumberOfDerivatives() const  {
  return 3*getNumberOfAtoms() + 9;
}

}
}

#endif
