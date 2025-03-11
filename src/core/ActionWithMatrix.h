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

class MatrixView;

class MatrixElementOutput {
public:
  View<double,helpers::dynamic_extent> values;
  View2D<double,helpers::dynamic_extent,helpers::dynamic_extent> derivs;
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
  void updateBookeepingArrays( MatrixView& outmat );
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithMatrix(const ActionOptions&);
  void performTask( const unsigned& task_index, MultiValue& myvals ) const override {
    plumed_error();
  }
/// Get the elements of the matrices into the output values
  void transferStashToValues( const std::vector<double>& stash ) override ;
/// Get the forces from the output values and transfer them to the stash
  void transferForcesToStash( std::vector<double>& stash ) const override ;
};

}
#endif
