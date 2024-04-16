/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#ifndef __PLUMED_dimred_SMACOF_h
#define __PLUMED_dimred_SMACOF_h

#include <vector>
#include "core/Value.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace dimred {

class SMACOF {
private:
  Matrix<double> Distances, Weights;
  double calculateSigma( const Matrix<double>& InitialZ, Matrix<double>& dists );
public:
  explicit SMACOF( const Value* mysquaredists );
  void optimize( const double& tol, const unsigned& maxloops, std::vector<double>& proj);
  double getDistance( const unsigned& i, const unsigned& j ) const ;
  void setWeight( const unsigned& i, const unsigned& j, const double& ww );
};

inline
double SMACOF::getDistance( const unsigned& i, const unsigned& j ) const {
  return Distances( i, j );
}

inline
void SMACOF::setWeight( const unsigned& i, const unsigned& j, const double& ww ) {
  Weights(i,j) = Weights(j,i ) = ww;
}

}
}
#endif
