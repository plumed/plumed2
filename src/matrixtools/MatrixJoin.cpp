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
#include "core/ActionRegister.h"

namespace PLMD {
namespace matrixtools {

class MatrixJoin : public ActionWithInputMatrices {
private:
  std::vector<unsigned> row_starts;
  std::vector<unsigned> col_starts;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit MatrixJoin(const ActionOptions&);
/// Do the calculation
  void completeMatrixOperations() override;
///
  void apply();
};

PLUMED_REGISTER_ACTION(MatrixJoin,"COMBINE_MATRICES")

void MatrixJoin::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrices::registerKeywords( keys ); keys.remove("ARG");
  keys.add("numbered","MATRIX","specify the matrices that you wish to join together into a single matrix"); keys.reset_style("MATRIX","compulsory");
}

MatrixJoin::MatrixJoin(const ActionOptions& ao):
  Action(ao),
  ActionWithInputMatrices(ao)
{
  plumed_assert( getNumberOfArguments()==0 );
  unsigned nrows=0, ncols=0; std::vector<Value*> arglist;
  for(unsigned i=1;; i++) {
    unsigned nt_cols=0; unsigned size_b4 = arglist.size();
    for(unsigned j=1;; j++) {
      if( j==10 ) error("cannot combine more than 10 matrices");
      std::vector<Value*> argn; parseArgumentList("MATRIX", i*10+j, argn);
      if( argn.size()==0 ) break;
      if( argn.size()>1 ) error("should only be one argument to each matrix keyword");
      if( argn[0]->getRank()!=2 ) error("input arguments for this action should be matrices");
      arglist.push_back( argn[0] ); nt_cols++;
      log.printf("  %d %d component of composed matrix is %s \n", i, j, argn[0]->getName().c_str() );
    }
    if( arglist.size()==size_b4 ) break;
    if( i==1 ) ncols=nt_cols;
    else if( nt_cols!=ncols ) error("should be joining same number of matrices in each row");
    nrows++;
  }

  std::vector<unsigned> shape(2); shape[0]=0; unsigned k=0;
  row_starts.resize( arglist.size() ); col_starts.resize( arglist.size() );
  for(unsigned i=0; i<nrows; ++i) {
    unsigned cstart = 0, nr = arglist[k]->getShape()[1];
    for(unsigned j=0; j<ncols; ++j) {
      if( arglist[k]->getShape()[1]!=nr ) error("mismatched matrix sizes");
      row_starts[k] = shape[0]; col_starts[k] = cstart;
      cstart += arglist[k]->getShape()[1]; k++;
    }
    if( i==0 ) shape[1]=cstart;
    else if( cstart!=shape[1] ) error("mismatched matrix sizes");
    shape[0] += arglist[k-1]->getShape()[0];
  }
  // Now request the arguments to make sure we store things we need
  requestArguments(arglist, false ); addValue( shape ); 
}

void MatrixJoin::completeMatrixOperations() {
  // Retrieve the matrix from input
  unsigned ncols = getPntrToOutput(0)->getShape()[1];
  for(unsigned k=0; k<getNumberOfArguments(); ++k) {
    Value* argn = getPntrToArgument(k); unsigned nn=0;
    // PLUMED cannot deal with sparsity when calculating these things
    plumed_assert( argn->getShape()[1]==(argn->getPntrToAction())->getNumberOfColumns() );
    for(unsigned i=0; i<argn->getShape()[0]; ++i) {
      for(unsigned j=0; j<argn->getShape()[1]; ++j) { getPntrToOutput(0)->set( ncols*(row_starts[k]+i)+col_starts[k]+j, argn->get(nn) ); nn++; }
    }
  }
}

void MatrixJoin::apply() {
  if( doNotCalculateDerivatives() ) return;

  if( getPntrToOutput(0)->forcesWereAdded() ) {
    unsigned ncols = getPntrToOutput(0)->getShape()[1];
    for(unsigned k=0; k<getNumberOfArguments(); ++k) {
      Value* argn = getPntrToArgument(k); unsigned nn=0;
      for(unsigned i=0; i<argn->getShape()[0]; ++i) {
        for(unsigned j=0; j<argn->getShape()[1]; ++j) { argn->addForce( nn, getPntrToOutput(0)->getForce(ncols*(row_starts[k]+i)+col_starts[k]+j) ); nn++; }
      }
    }
  }

}

}
}
