/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#ifndef __PLUMED_tools_LatticeReduction_h
#define __PLUMED_tools_LatticeReduction_h

#include "Vector.h"
#include "Tensor.h"

namespace PLMD {

/**
Class implementing algorithms for lattice reduction.

This class implements algorithms described in
Igor Semaev, A 3-Dimensional Lattice Reduction Algorithm, CaLC 2001, LNCS 2146, pp. 181–193, 2001.
It just collect static methods in a separate namespace.
*/
class LatticeReduction {
/// Gaussian reduction
  static void reduce(Vector&a,Vector&b);
/// Obtain three reduce-2 vectors (Algorithm 1 in the paper), equivalent to reduce2(Tensor&t)
  static void reduce2(Vector&a,Vector&b,Vector&c);
/// Check if two vectors are reduced
  static bool isReduced(const Vector&a,const Vector&b);
/// Check if three vectors are reduced
  static bool isReduced(const Vector&a,const Vector&b,const Vector &c);
/// Check if three vectors are reduced-2
  static bool isReduced2(const Vector&a,const Vector&b,const Vector &c);
/// Obtain three reduce-2 vectors (Algorithm 1 in the paper), equivalent to reduce2(Vector&a,Vector&b,Vector&c)
  static void reduce2(Tensor&t);
/// Sort three vectors by modulo
  static void sort(Vector v[3]);
public:
/// Reduce a basis in place, maps to reduceFast()
  static void reduce(Tensor&t);
/// Reduce a basis in place using the slow algorithm (Algorithm 2 in the paper)
  static void reduceSlow(Tensor&t);
/// Reduce a basis in place using the fast algorithm (Algorithm 3 in the paper)
  static void reduceFast(Tensor&t);
/// Check if a basis is reduced
  static bool isReduced(const Tensor&t);
};

}

#endif

