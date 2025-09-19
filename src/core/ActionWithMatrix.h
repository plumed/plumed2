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
#ifndef __PLUMED_core_ActionWithMatrix_h
#define __PLUMED_core_ActionWithMatrix_h

#include "ActionWithVector.h"

namespace PLMD {

class ActionWithMatrix;
class RequiredMatrixElements {
  std::vector<std::size_t> bookeeping;
  std::size_t const* bookeeping_data;
public:
  friend class ActionWithMatrix;
  std::size_t ncols;
  void update() {
    bookeeping_data = bookeeping.data();
  }
  void toACCDevice() const {
    //this assumes that update() has already been called
#pragma acc enter data copyin(this[0:1], bookeeping_data[0:bookeeping.size()])
  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(bookeeping_data[0:bookeeping.size()],this[0:1])
  }
  std::size_t operator[]( std::size_t i ) const {
    return bookeeping_data[i];
  }

};

class MatrixElementOutput {
public:
  View<double> values;
  View2D<double> derivs;
  MatrixElementOutput( std::size_t nv, std::size_t nd, double* v, double* d ) :
    values( v, nv ),
    derivs( d, nv, nd ) {
  }
};

class ActionWithMatrix : public ActionWithVector {
protected:
/// Some actions have a flag that sets this to true.  The elements on the diagonal of the resulting matrix are then set to zero
  bool diagzero;
/// Update all the arrays for doing bookeeping
  void updateBookeepingArrays( RequiredMatrixElements& outmat );
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithMatrix(const ActionOptions&);
/// Get the elements of the matrices into the output values
  void transferStashToValues( const std::vector<double>& stash ) override ;
/// Get the elements of the matrices into the output values
  void transferStashToValues( const std::vector<float>& stash ) override ;
/// Get the forces from the output values and transfer them to the stash
  void transferForcesToStash( std::vector<double>& stash ) const override ;
/// Get the forces from the output values and transfer them to the stash
  void transferForcesToStash( std::vector<float>& stash ) const override ;
};

}
#endif
