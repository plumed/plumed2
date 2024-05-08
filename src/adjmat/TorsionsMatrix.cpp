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
#include "core/ActionWithMatrix.h"
#include "core/ActionRegister.h"
#include "tools/Torsion.h"

//+PLUMEDOC MCOLVAR TORSIONS_MATRIX
/*
Calculate the matrix of torsions between two vectors of molecules

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class TorsionsMatrix : public ActionWithMatrix {
private:
  unsigned nderivatives;
  bool stored_matrix1, stored_matrix2;
public:
  static void registerKeywords( Keywords& keys );
  explicit TorsionsMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  unsigned getNumberOfColumns() const override { return getConstPntrToComponent(0)->getShape()[1]; }
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override;
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override ;
};

PLUMED_REGISTER_ACTION(TorsionsMatrix,"TORSIONS_MATRIX")

void TorsionsMatrix::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys); keys.use("ARG");
  keys.add("atoms","POSITIONS1","the positions to use for the molecules specified using the first argument");
  keys.add("atoms","POSITIONS2","the positions to use for the molecules specified using the second argument");
  keys.setValueDescription("the matrix of torsions between the two vectors of input directors");
}

TorsionsMatrix::TorsionsMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao)
{
  if( getNumberOfArguments()!=2 ) error("should be two arguments to this action, a matrix and a vector");
  if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) error("first argument to this action should be a matrix");
  if( getPntrToArgument(1)->getRank()!=2 || getPntrToArgument(1)->hasDerivatives() ) error("second argument to this action should be a matrix");
  if( getPntrToArgument(0)->getShape()[1]!=3 || getPntrToArgument(1)->getShape()[0]!=3 ) error("number of columns in first matrix and number of rows in second matrix should equal 3");

  std::vector<AtomNumber> atoms_a; parseAtomList("POSITIONS1", atoms_a );
  if( atoms_a.size()!=getPntrToArgument(0)->getShape()[0] ) error("mismatch between number of atoms specified using POSITIONS1 and number of arguments in vector input");
  log.printf("  using positions of these atoms for vectors in first matrix \n");
  for(unsigned int i=0; i<atoms_a.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", atoms_a[i].serial());
  }
  log.printf("\n"); std::vector<AtomNumber> atoms_b; parseAtomList("POSITIONS2", atoms_b );
  if( atoms_b.size()!=getPntrToArgument(1)->getShape()[1] ) error("mismatch between number of atoms specified using POSITIONS2 and number of arguments in vector input");
  log.printf("  using positions of these atoms for vectors in second matrix \n");
  for(unsigned i=0; i<atoms_b.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", atoms_b[i].serial()); atoms_a.push_back( atoms_b[i] );
  }
  log.printf("\n"); requestAtoms( atoms_a, false );

  std::vector<unsigned> shape(2); shape[0]=getPntrToArgument(0)->getShape()[0]; shape[1]=getPntrToArgument(1)->getShape()[1];
  addValue( shape ); setPeriodic("-pi","pi"); nderivatives = buildArgumentStore(0) + 3*getNumberOfAtoms() + 9;
  std::string headstr=getFirstActionInChain()->getLabel();
  stored_matrix1 = getPntrToArgument(0)->ignoreStoredValue( headstr );
  stored_matrix2 = getPntrToArgument(1)->ignoreStoredValue( headstr );
}

unsigned TorsionsMatrix::getNumberOfDerivatives() {
  return nderivatives;
}

void TorsionsMatrix::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  unsigned start_n = getPntrToArgument(0)->getShape()[0], size_v = getPntrToArgument(1)->getShape()[1];
  if( indices.size()!=size_v+1 ) indices.resize( size_v+1 );
  for(unsigned i=0; i<size_v; ++i) indices[i+1] = start_n + i;
  myvals.setSplitIndex( size_v + 1 );
}

void TorsionsMatrix::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned ostrn = getConstPntrToComponent(0)->getPositionInStream(), ind2=index2;
  if( index2>=getPntrToArgument(0)->getShape()[0] ) ind2 = index2 - getPntrToArgument(0)->getShape()[0];

  Vector v1, v2, dv1, dv2, dconn;
  // Compute the distance connecting the two centers
  Vector conn=pbcDistance( getPosition(index1), getPosition(index2) );
  if( conn.modulo2()<epsilon ) return;

  // Get the two vectors
  for(unsigned i=0; i<3; ++i) {
    v1[i] = getElementOfMatrixArgument( 0, index1, i, myvals );
    v2[i] = getElementOfMatrixArgument( 1, i, ind2, myvals );
  }
  // Evaluate angle
  Torsion t; double angle = t.compute( v1, conn, v2, dv1, dconn, dv2 );
  myvals.addValue( ostrn, angle );

  if( doNotCalculateDerivatives() ) return;

  // Add the derivatives on the matrices
  for(unsigned i=0; i<3; ++i) {
    addDerivativeOnMatrixArgument( stored_matrix1, 0, 0, index1, i, dv1[i], myvals );
    addDerivativeOnMatrixArgument( stored_matrix2, 0, 1, i, ind2, dv2[i], myvals );
  }
  // And derivatives on positions
  unsigned narg_derivatives = getPntrToArgument(0)->getNumberOfValues() + getPntrToArgument(1)->getNumberOfValues();
  for(unsigned i=0; i<3; ++i) {
    myvals.addDerivative( ostrn, narg_derivatives + 3*index1+i, -dconn[i] ); myvals.addDerivative( ostrn, narg_derivatives + 3*index2+i, dconn[i] );
    myvals.updateIndex( ostrn, narg_derivatives + 3*index1+i ); myvals.updateIndex( ostrn, narg_derivatives + 3*index2+i );
  }
  //And virial
  Tensor vir( -extProduct( conn, dconn ) ); unsigned virbase = narg_derivatives + 3*getNumberOfAtoms();
  for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j ) { myvals.addDerivative( ostrn, virbase+3*i+j, vir(i,j) ); myvals.updateIndex( ostrn, virbase+3*i+j ); }
}

void TorsionsMatrix::runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() || !matrixChainContinues() ) return ;

  unsigned mat1s = 3*ival, ss = getPntrToArgument(1)->getShape()[1];
  unsigned nmat = getConstPntrToComponent(0)->getPositionInMatrixStash(), nmat_ind = myvals.getNumberOfMatrixRowDerivatives( nmat );
  unsigned narg_derivatives = getPntrToArgument(0)->getNumberOfValues() + getPntrToArgument(1)->getNumberOfValues();
  std::vector<unsigned>& matrix_indices( myvals.getMatrixRowDerivativeIndices( nmat ) ); unsigned ntwo_atoms = myvals.getSplitIndex();
  for(unsigned j=0; j<3; ++j) {
    matrix_indices[nmat_ind] = mat1s + j; nmat_ind++;
    matrix_indices[nmat_ind] = narg_derivatives + mat1s + j; nmat_ind++;
    for(unsigned i=1; i<ntwo_atoms; ++i) {
      unsigned ind2 = indices[i]; if( ind2>=getPntrToArgument(0)->getShape()[0] ) ind2 = indices[i] - getPntrToArgument(0)->getShape()[0];
      matrix_indices[nmat_ind] = arg_deriv_starts[1] + j*ss + ind2; nmat_ind++;
      matrix_indices[nmat_ind] = narg_derivatives + 3*indices[i] + j; nmat_ind++;
    }
  }
  unsigned base = narg_derivatives + 3*getNumberOfAtoms(); for(unsigned j=0; j<9; ++j) { matrix_indices[nmat_ind] = base + j; nmat_ind++; }
  myvals.setNumberOfMatrixRowDerivatives( nmat, nmat_ind );
}

}
}
