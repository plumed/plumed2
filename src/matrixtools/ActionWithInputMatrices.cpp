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
  for(unsigned i=0; i<getNumberOfArguments();++i) {
      if( getPntrToArgument(i)->getRank()!=2 ) error("input argument for this action should be a matrix");
  }
  // Now request the arguments to make sure we store things we need
  std::vector<Value*> args( getArguments() ); arg_ends.push_back(0); arg_ends.push_back(args.size()); 
  requestArguments(args, false ); 
}

void ActionWithInputMatrices::addValue( const std::vector<unsigned>& shape ) {
  ActionWithValue::addValue( shape ); setNotPeriodic(); getPntrToOutput(0)->alwaysStoreValues();
  bool istimeseries=false;
  for(unsigned i=0; i<getNumberOfArguments();++i) {
      if( getPntrToArgument(0)->isTimeSeries() ) { istimeseries=true; break; }
  }

  if( istimeseries  ) {
      for(unsigned i=0; i<getNumberOfArguments();++i) {
          if( !getPntrToArgument(0)->isTimeSeries() ) error( "on argument is time series but " + getPntrToArgument(i)->getName() + " is not a time series");
      }
      getPntrToOutput(0)->makeTimeSeries();
  }
}

unsigned ActionWithInputMatrices::getNumberOfDerivatives() const {
  return 0;
}

unsigned ActionWithInputMatrices::getNumberOfColumns() const {
  return getPntrToOutput(0)->getShape()[1];
}

void ActionWithInputMatrices::retrieveFullMatrix( const unsigned& imat, Matrix<double>& mymatrix ) const {
  unsigned k = 0;
  for(unsigned i=0; i<mymatrix.nrows(); ++i) {
    for(unsigned j=0; j<mymatrix.ncols(); ++j) { mymatrix(i,j) = getPntrToArgument(imat)->get( k ); k++; }
  }
}

void ActionWithInputMatrices::calculate() {
  if( skipCalculate() ) return;
  completeMatrixOperations();
}

void ActionWithInputMatrices::update() {
  if( skipUpdate() ) return;
  completeMatrixOperations();
}

void ActionWithInputMatrices::runFinalJobs() {
  if( skipUpdate() ) return;
  if( getName()=="VORONOI" ) {
      std::vector<unsigned> shape( getPntrToArgument(0)->getShape() );
      for(unsigned i=0; i<shape[1]; ++i) addTaskToList( i );
      getPntrToOutput(0)->setShape( shape );
  } else resizeForFinalTasks(); 
  completeMatrixOperations();
}

}
}
