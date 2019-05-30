/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#ifndef __PLUMED_tools_Pbc_h
#define __PLUMED_tools_Pbc_h

#include "Vector.h"
#include "Tensor.h"
#include <vector>
#include <cstddef>

namespace PLMD {

/*
Tool to deal with periodic boundary conditions.

This class is useful to apply periodic boundary conditions on interatomic
distances. It stores privately information about reduced lattice vectors
*/
class Pbc {
/// Type of box
  enum {unset,orthorombic,generic} type;
/// Box
  Tensor box;
/// Inverse box
  Tensor invBox;
/// Reduced box.
/// This is a set of lattice vectors generating the same lattice
/// but "minimally skewed". Useful to optimize image search.
  Tensor reduced;
/// Inverse of the reduced box
  Tensor invReduced;
/// List of shifts that should be attempted.
/// Depending on the sign of the scaled coordinates representing
/// a distance vector, a different set of shifts must be tried.
  std::vector<Vector> shifts[2][2][2];
/// Alternative representation for orthorombic cells.
/// Not really used, but could be used to optimize search in
/// orthorombic cells.
  Vector diag,hdiag,mdiag;
/// Build list of shifts.
/// This is expensive, and must be called only when box is
/// reset. It allows building a minimal set of shifts
/// depending on the sign of the scaled coordinates representing
/// a distance vector.
  void buildShifts(std::vector<Vector> shifts[2][2][2])const;
public:
/// Constructor
  Pbc();
/// Compute modulo of (v2-v1), using or not pbc depending on bool pbc.
  double distance( const bool pbc, const Vector& v1, const Vector& v2 ) const;
/// Computes v2-v1, using minimal image convention
  Vector distance(const Vector& v1,const Vector& v2)const;
/// version of distance which also returns the number
/// of attempted shifts
  Vector distance(const Vector&,const Vector&,int*nshifts)const;
/// Apply PBC to a set of positions or distance vectors
  void apply(std::vector<Vector>&dlist, unsigned max_index=0) const;
/// Set the lattice vectors.
/// b[i][j] is the j-th component of the i-th vector
  void setBox(const Tensor&b);
/// Returns the box
  const Tensor& getBox()const;
/// Returns the inverse matrix of box.
/// Thus: pbc.getInvBox() == inverse(pbc.getBox()).
  const Tensor& getInvBox()const;
/// Transform a vector in real space to a vector in scaled coordinates.
/// Thus:pbc.realToScaled(v) == matmul(transpose(inverse(pbc.getBox(),v)));
  Vector realToScaled(const Vector&)const;
/// Transform a vector in scaled coordinates to a vector in real space.
/// Thus:pbc.scaledToRead(v) == matmul(transpose(pbc.getBox()),v);
  Vector scaledToReal(const Vector&)const;
/// Returns true if the box vectors are orthogonal
  bool isOrthorombic()const;
/// Full search (for testing).
/// Perform a full search on vector
  void fullSearch(Vector&)const;
/// Returns true if box is set and non zero
  bool isSet()const;
};

inline
Vector Pbc::distance(const Vector& v1,const Vector& v2)const {
  return distance(v1,v2,NULL);
}

inline
bool Pbc::isSet()const {
  return type!=unset;
}

}

#endif
