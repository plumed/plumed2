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
#include "ClassicalScaling.h"
#include "reference/PointWiseMapping.h"

namespace PLMD {
namespace analysis {

void ClassicalScaling::run( PointWiseMapping* mymap ) {
  // Retrieve the distances from the dimensionality reduction object
  double half=(-0.5); Matrix<double> distances( half*mymap->modifyDmat() );

  // Apply centering transtion
  unsigned n=distances.nrows(); double sum;
  // First HM
  for(unsigned i=0; i<n; ++i) {
    sum=0; for(unsigned j=0; j<n; ++j) sum+=distances(i,j);
    for(unsigned j=0; j<n; ++j) distances(i,j) -= sum/n;
  }
  // Now (HM)H
  for(unsigned i=0; i<n; ++i) {
    sum=0; for(unsigned j=0; j<n; ++j) sum+=distances(j,i);
    for(unsigned j=0; j<n; ++j) distances(j,i) -= sum/n;
  }

  // Diagonalize matrix
  std::vector<double> eigval(n); Matrix<double> eigvec(n,n);
  diagMat( distances, eigval, eigvec );

  // Pass final projections to map object
  for(unsigned i=0; i<n; ++i) {
    for(unsigned j=0; j<mymap->getNumberOfProperties(); ++j) mymap->setProjectionCoordinate( i, j, sqrt(eigval[n-1-j])*eigvec(n-1-j,i) );
  }
}

}
}
