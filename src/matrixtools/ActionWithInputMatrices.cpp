/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "ActionWithInputMatrices.h"

namespace PLMD {
namespace matrixtools {

void ActionWithInputMatrices::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithArguments::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); keys.use("ARG"); 
}

ActionWithInputMatrices::ActionWithInputMatrices(const ActionOptions& ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao)
{
  if( getName()!="TRANSPOSE" && getName()!="SELECT_COMPONENTS" && getName()!="CONCATENATE" ) {
      for(unsigned i=0; i<getNumberOfArguments();++i) {
          if( getPntrToArgument(i)->getRank()!=2 ) error("input argument for this action should be a matrix");
      }
  }
  // Now request the arguments to make sure we store things we need
  std::vector<Value*> args( getArguments() ); 
  requestArguments(args, false ); 
}

void ActionWithInputMatrices::addValue( const std::vector<unsigned>& shape ) {
  ActionWithValue::addValue( shape ); 
  if( getName()=="TRANSPOSE" && getPntrToArgument(0)->isPeriodic() ) { 
      std::string smin, smax; getPntrToArgument(0)->getDomain( smin, smax ); setPeriodic( smin, smax );
  } else setNotPeriodic();
  getPntrToOutput(0)->alwaysStoreValues();
}

unsigned ActionWithInputMatrices::getNumberOfDerivatives() const {
  return 0;
}

unsigned ActionWithInputMatrices::getNumberOfColumns() const {
  return getPntrToOutput(0)->getShape()[0];
}

void ActionWithInputMatrices::retrieveEdgeList( const unsigned& imat, unsigned& nedge ) {
  nedge=0; Value* mat = getPntrToArgument(imat); 
  unsigned nrows = mat->getShape()[0], ncols = mat->getNumberOfColumns();
  // Check we have enough space to store the edge list
  if( pairs.size()<nrows*ncols ) { vals.resize( nrows*ncols ); pairs.resize( nrows*ncols ); }

  bool symmetric = mat->isSymmetric();
  for(unsigned i=0; i<nrows; ++i) {
      unsigned ncol = mat->getRowLength(i);
      for(unsigned j=0; j<ncol; ++j) {
          if( fabs(mat->get(i*ncols+j,false))<epsilon ) continue;
          if( symmetric && mat->getRowIndex(i,j)>i ) continue;
          pairs[nedge].first = i; pairs[nedge].second = mat->getRowIndex(i,j); 
          vals[nedge] = mat->get(i*ncols+j,false); nedge++;
      }
  }
}

void ActionWithInputMatrices::retrieveFullMatrix( const unsigned& imat, Matrix<double>& mymatrix ) {
  unsigned nedge; retrieveEdgeList( imat, nedge  ); mymatrix=0; 
  bool symmetric=getPntrToArgument(imat)->isSymmetric();
  for(unsigned i=0; i<nedge; ++i ) {
      mymatrix( pairs[i].first, pairs[i].second ) = vals[i];  
      if( symmetric ) mymatrix( pairs[i].second, pairs[i].first ) = vals[i];
  } 
}

void ActionWithInputMatrices::calculate() {
  if( skipCalculate() ) return;
  completeMatrixOperations();
}

void ActionWithInputMatrices::applyForceOnMatrix( const unsigned& imat ) {
  Value* mat = getPntrToArgument(imat);

  unsigned ncols=mat->getNumberOfColumns();
  for(unsigned i=0; i<mat->getShape()[0]; ++i) {
      unsigned ncol = mat->getRowLength(i);
      for(unsigned j=0; j<ncol; ++j) mat->addForce( i*ncols+j, getForceOnMatrixElement( imat, i, mat->getRowIndex(i,j) ), false );
  }
}

void ActionWithInputMatrices::update() {
  if( skipUpdate() ) return;
  completeMatrixOperations();
}

void ActionWithInputMatrices::runFinalJobs() {
  if( skipUpdate() ) return;
  if( getName()=="VORONOI" ) {
      std::vector<unsigned> shape( getPntrToArgument(0)->getShape() );
      getPntrToOutput(0)->setShape( shape ); getPntrToOutput(0)->setNumberOfTasks( shape[1] );
  } else if( getName()!="SELECT_COMPONENTS" && getName()!="DIAGONALIZE" ) getPntrToOutput(0)->setShape( getValueShapeFromArguments() ); 
  completeMatrixOperations();
}

}
}
