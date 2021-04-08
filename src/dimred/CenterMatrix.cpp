/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "matrixtools/ActionWithInputMatrices.h"
#include "core/ActionRegister.h"

//+PLUMEDOC DIMRED CENTER_MATRIX
/*
Create a low-dimensional projection of a trajectory using the classical multidimensional
 scaling algorithm.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class CenterMatrix : public matrixtools::ActionWithInputMatrices {
private:
  bool wasupdated;
public:
  static void registerKeywords( Keywords& keys );
  explicit CenterMatrix( const ActionOptions& ao );
  void completeMatrixOperations() override;
  void apply(){}
};

PLUMED_REGISTER_ACTION(CenterMatrix,"CENTER_MATRIX")

void CenterMatrix::registerKeywords( Keywords& keys ) {
  matrixtools::ActionWithInputMatrices::registerKeywords( keys ); 
}

CenterMatrix::CenterMatrix( const ActionOptions& ao):
  Action(ao),
  matrixtools::ActionWithInputMatrices(ao),
  wasupdated(false)
{
  if( getNumberOfArguments()!=1 ) error("should only input one argument to CENTER_MATRIX");
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) error("input to CENTER_MATRIX should be a square matrix");

  // Re-request the arguments as we do in diagonalize for similar reasons
  std::vector<unsigned> shape(2); shape[0]=shape[1] = getPntrToArgument(0)->getShape()[0];
  addValue(shape); setNotPeriodic(); getPntrToOutput(0)->alwaysStoreValues();
  if( getPntrToArgument(0)->isTimeSeries() ) getPntrToOutput(0)->makeTimeSeries();
}

void CenterMatrix::completeMatrixOperations() {
  plumed_dbg_assert( !actionInChain() ); Value* arg0=getPntrToArgument(0); Value* val0=getPntrToOutput(0);

  // Apply centering transtion
  unsigned n=arg0->getShape()[0]; 
  if( n==0 ) return ;

  // First HM
  double sum; wasupdated=true;
  for(unsigned i=0; i<n; ++i) {
    sum=0; for(unsigned j=0; j<n; ++j) sum+=-0.5*arg0->get( i*n+j );             
    sum /= n; for(unsigned j=0; j<n; ++j) val0->set( i*n+j, -0.5*arg0->get( i*n+j ) - sum );  
  }
  // Now (HM)H
  for(unsigned i=0; i<n; ++i) {
    sum=0; for(unsigned j=0; j<n; ++j) sum+=val0->get( j*n+i );    
    sum /= n; for(unsigned j=0; j<n; ++j) val0->set( j*n+i, val0->get( j*n+i ) - sum ); 
  }
}

}
}
