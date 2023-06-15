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
#include "MatrixOperationBase.h"

namespace PLMD {
namespace adjmat {

void MatrixOperationBase::retrieveEdgeList( const Value* mat, unsigned& nedge, std::vector<std::pair<unsigned,unsigned> >& active, std::vector<double>& elems ) {
  nedge=0; plumed_dbg_assert( mat->getRank()==2 && !mat->hasDerivatives() );
  unsigned nrows = mat->getShape()[0], ncols = mat->getNumberOfColumns();
  // Check we have enough space to store the edge list
  if( elems.size()<nrows*ncols ) { elems.resize( nrows*ncols ); active.resize( nrows*ncols ); }

  bool symmetric = mat->isSymmetric();
  for(unsigned i=0; i<nrows; ++i) {
      unsigned ncol = mat->getRowLength(i);
      for(unsigned j=0; j<ncol; ++j) {
          if( fabs(mat->get(i*ncols+j,false))<epsilon ) continue;
          if( symmetric && mat->getRowIndex(i,j)>i ) continue;
          active[nedge].first = i; active[nedge].second = mat->getRowIndex(i,j);
          elems[nedge] = mat->get(i*ncols+j,false); nedge++;
      }
  }
}

void MatrixOperationBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithArguments::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  keys.use("ARG"); keys.remove("NUMERICAL_DERIVATIVES");
}

MatrixOperationBase::MatrixOperationBase(const ActionOptions&ao):
Action(ao),
ActionWithArguments(ao),
ActionWithValue(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument to this action");
  if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) error("input to this argument should be a matrix");
  getPntrToArgument(0)->buildDataStore();
}

void MatrixOperationBase::retrieveFullMatrix( Matrix<double>& mymatrix ) {
  if( mymatrix.nrows()!=getPntrToArgument(0)->getShape()[0] || mymatrix.nrows()!=getPntrToArgument(0)->getShape()[1] ) {
      mymatrix.resize( getPntrToArgument(0)->getShape()[0], getPntrToArgument(0)->getShape()[1] );
  }
  unsigned nedge; retrieveEdgeList( getPntrToArgument(0), nedge, pairs, vals  ); mymatrix=0;
  bool symmetric = getPntrToArgument(0)->isSymmetric();
  for(unsigned i=0; i<nedge; ++i ) {
      mymatrix( pairs[i].first, pairs[i].second ) = vals[i];
      if( symmetric ) mymatrix( pairs[i].second, pairs[i].first ) = vals[i];
  }
}


void MatrixOperationBase::apply() {
  if( doNotCalculateDerivatives() ) return;

  bool forces=false;
  for(int i=0; i<getNumberOfComponents();++i) {
      if( getPntrToComponent(i)->forcesWereAdded() ) { forces=true; break; }
  }
  if( !forces ) return;

  Value* mat = getPntrToArgument(0);
  unsigned ncols=mat->getNumberOfColumns();
  for(unsigned i=0; i<mat->getShape()[0]; ++i) {
      unsigned ncol = mat->getRowLength(i);
      for(unsigned j=0; j<ncol; ++j) mat->addForce( i*ncols+j, getForceOnMatrixElement( i, mat->getRowIndex(i,j) ), false );
  }
}


}
}
