/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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
#include "MatrixView.h"

namespace PLMD {

void MatrixView::setup( std::size_t s, const Value* myval ) {
  start = s;
  shape.resize(2);
  if( myval->getRank()==1 ) {
    shape[0] = myval->getShape()[0];
    shape[1] = ncols = 1;
    bookeeping.resize( 2*shape[0] );
    for(unsigned i=0; i<shape[0]; ++i) {
      bookeeping[2*i] = 1;
      bookeeping[2*i+1] = 0;
    }
  } else if( myval->getRank()==2 ) {
    shape[0] = myval->getShape()[0];
    shape[1] = myval->getShape()[1];
    ncols = myval->getNumberOfColumns();
    bookeeping.resize( myval->matrix_bookeeping.size() );
    for(unsigned i=0; i<bookeeping.size(); ++i) {
      bookeeping[i] = myval->matrix_bookeeping[i];
    }
  }
}

double MatrixView::getElement( std::size_t irow, std::size_t jcol, const MatrixView& mat, std::vector<double>& data ) {
  if( mat.shape[1]==mat.ncols ) {
    plumed_assert( mat.start + irow*mat.ncols + jcol<data.size() );
    return data[mat.start + irow*mat.ncols + jcol];
  }

  for(unsigned i=0; i<mat.bookeeping[(1+mat.ncols)*irow]; ++i) {
    if( mat.bookeeping[(1+mat.ncols)*irow+1+i]==jcol ) {
      return data[mat.start + irow*mat.ncols+i];
    }
  }
  return 0.0;
}

bool MatrixView::hasElement( std::size_t irow, std::size_t jcol, const MatrixView& mat, std::size_t& ind ) {
  if( mat.shape[1]==mat.ncols ) {
    ind = mat.start + irow*mat.ncols + jcol;
    return true;
  }
  for(unsigned i=0; i<mat.bookeeping[(1+mat.ncols)*irow]; ++i) {
    if( mat.bookeeping[(1+mat.ncols)*irow+1+i]==jcol ) {
      ind = mat.start + irow*mat.ncols+i;
      return true;
    }
  }
  return false;
}

}
