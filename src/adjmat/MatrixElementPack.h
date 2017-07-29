/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#ifndef __PLUMED_adjmat_MatrixElementPack_h
#define __PLUMED_adjmat_MatrixElementPack_h

#include <vector>
#include "tools/MultiValue.h"
#include "tools/Vector.h"
#include "tools/Tensor.h"

namespace PLMD {
namespace adjmat {

class MatrixElementPack {
private:
/// Copy of the values that we are adding to
  MultiValue& myvals;
///
  unsigned natoms_in_base;
///
  bool noderiv;
/// Index for weight of this connection
  unsigned w_index;
/// Indices for x, y and z components of vector
  bool components;
  unsigned x_index, y_index, z_index;
///
  std::vector<unsigned> indices;
///
  std::vector<Vector> positions;
public:
  MatrixElementPack( MultiValue& vals, const bool comp, ActionWithValue const * myaction );
  void setIndex( const unsigned& indno, const unsigned& ind );
  void setPosition( const unsigned& indno, const Vector& pos );
  Vector getPosition( const unsigned& indno ) const ;
  void addAtomDerivatives( const unsigned& indno, const Vector& der );
  void addBoxDerivatives( const Tensor& vir );
  void clear();
};

inline
void MatrixElementPack::setIndex( const unsigned& indno, const unsigned& ind ){
  plumed_dbg_assert( indno<indices.size() ); indices[indno]=ind;
}

inline
void MatrixElementPack::setPosition( const unsigned& indno, const Vector& pos ){
  plumed_dbg_assert( natoms_in_base>0 && indno<indices.size() ); 
  positions[indno]=pos;
}

inline
Vector MatrixElementPack::getPosition( const unsigned& indno ) const {
  plumed_dbg_assert( natoms_in_base>0 && indno<indices.size() ); 
  return positions[indno];
}

inline
void MatrixElementPack::addAtomDerivatives( const unsigned& indno, const Vector& der ){
  if( noderiv ) return;
  plumed_dbg_assert( natoms_in_base>0 && indno<indices.size() );
  myvals.addDerivative( w_index, 3*indices[indno]+0, der[0] );
  myvals.addDerivative( w_index, 3*indices[indno]+1, der[1] );
  myvals.addDerivative( w_index, 3*indices[indno]+2, der[2] );
}

inline
void MatrixElementPack::addBoxDerivatives( const Tensor& vir ){
  if( noderiv ) return;
  plumed_dbg_assert( natoms_in_base>0 );
  myvals.addDerivative( w_index, 3*natoms_in_base+0, vir(0,0) );
  myvals.addDerivative( w_index, 3*natoms_in_base+1, vir(0,1) );
  myvals.addDerivative( w_index, 3*natoms_in_base+2, vir(0,2) );
  myvals.addDerivative( w_index, 3*natoms_in_base+3, vir(1,0) );
  myvals.addDerivative( w_index, 3*natoms_in_base+4, vir(1,1) );
  myvals.addDerivative( w_index, 3*natoms_in_base+5, vir(1,2) );
  myvals.addDerivative( w_index, 3*natoms_in_base+6, vir(2,0) );
  myvals.addDerivative( w_index, 3*natoms_in_base+7, vir(2,1) );
  myvals.addDerivative( w_index, 3*natoms_in_base+8, vir(2,2) );
}


}
}
#endif
