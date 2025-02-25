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
public:
  static void registerKeywords( Keywords& keys );
  explicit TorsionsMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  unsigned getNumberOfColumns() const override ;
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override;
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override ;
};

PLUMED_REGISTER_ACTION(TorsionsMatrix,"TORSIONS_MATRIX")

void TorsionsMatrix::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys); 
  keys.addInputKeyword("optional","MASK","matrix","a matrix that is used to used to determine which elements of this matrix to compute");
  keys.addInputKeyword("compulsory","ARG","matrix","an Nx3 and a 3xN matrix that contain the bond vectors that you would like to determine the torsion angles between");
  keys.add("atoms","POSITIONS1","the positions to use for the molecules specified using the first argument");
  keys.add("atoms","POSITIONS2","the positions to use for the molecules specified using the second argument");
  keys.setValueDescription("matrix","the matrix of torsions between the two vectors of input directors");
}

TorsionsMatrix::TorsionsMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao)
{
  int nm=getNumberOfMasks(); if( nm<0 ) nm = 0;
  if( getNumberOfArguments()-nm!=2 ) error("should be two arguments to this action, a matrix and a vector");
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
  addValue( shape ); setPeriodic("-pi","pi");

  if( nm>0 ) {
    unsigned iarg = getNumberOfArguments()-1;
    if( getPntrToArgument(iarg)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) error("argument passed to MASK keyword should be a matrix");
    if( getPntrToArgument(iarg)->getShape()[0]!=shape[0] || getPntrToArgument(iarg)->getShape()[1]!=shape[1] ) error("argument passed to MASK keyword has the wrong shape");
  }
}

unsigned TorsionsMatrix::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms() + 9 + getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(1)->getNumberOfStoredValues();
}

unsigned TorsionsMatrix::getNumberOfColumns() const {
  if( getNumberOfArguments()>2 ) return getPntrToArgument(2)->getNumberOfColumns();
  return getConstPntrToComponent(0)->getShape()[1];
}

void TorsionsMatrix::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  if( getNumberOfArguments()>2 ) {
    unsigned start_n = getPntrToArgument(0)->getShape()[0], size_v = getPntrToArgument(2)->getRowLength(task_index);
    if( indices.size()!=size_v+1 ) indices.resize( size_v+1 );
    for(unsigned i=0; i<size_v; ++i) indices[i+1] = start_n + getPntrToArgument(2)->getRowIndex(task_index, i);
    myvals.setSplitIndex( size_v + 1 );
  } else {
    unsigned start_n = getPntrToArgument(0)->getShape()[0], size_v = getPntrToArgument(1)->getShape()[1];
    if( indices.size()!=size_v+1 ) indices.resize( size_v+1 );
    for(unsigned i=0; i<size_v; ++i) indices[i+1] = start_n + i;
    myvals.setSplitIndex( size_v + 1 );
  }
}

void TorsionsMatrix::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned ind2=index2;
  if( index2>=getPntrToArgument(0)->getShape()[0] ) ind2 = index2 - getPntrToArgument(0)->getShape()[0];

  Vector v1, v2, dv1, dv2, dconn;
  // Compute the distance connecting the two centers
  Vector conn=pbcDistance( getPosition(index1), getPosition(index2) );
  if( conn.modulo2()<epsilon ) return;

  // Get the two vectors
  for(unsigned i=0; i<3; ++i) {
    v1[i] = getPntrToArgument(0)->get( index1*getPntrToArgument(0)->getNumberOfColumns() + i, false );
    v2[i] = getPntrToArgument(1)->get( i*getPntrToArgument(1)->getNumberOfColumns() + ind2, false );
  }
  // Evaluate angle
  Torsion t; double angle = t.compute( v1, conn, v2, dv1, dconn, dv2 );
  myvals.addValue( 0, angle );

  if( doNotCalculateDerivatives() ) return;

  // Add the derivatives on the matrices
  unsigned base1 = index1*getPntrToArgument(0)->getNumberOfColumns();
  unsigned ncols = getPntrToArgument(1)->getNumberOfColumns();
  unsigned base2 = getPntrToArgument(0)->getNumberOfStoredValues() + ind2;
  for(unsigned i=0; i<3; ++i) {
    myvals.addDerivative( 0, base1 + i, dv1[i] ); myvals.updateIndex( 0, base1 + i );
    myvals.addDerivative( 0, base2 + i*ncols, dv2[i] ); myvals.updateIndex( 0, base2 + i*ncols );
  }
  // And derivatives on positions
  unsigned narg_derivatives = getPntrToArgument(0)->getNumberOfValues() + getPntrToArgument(1)->getNumberOfValues();
  for(unsigned i=0; i<3; ++i) {
    myvals.addDerivative( 0, narg_derivatives + 3*index1+i, -dconn[i] ); myvals.addDerivative( 0, narg_derivatives + 3*index2+i, dconn[i] );
    myvals.updateIndex( 0, narg_derivatives + 3*index1+i ); myvals.updateIndex( 0, narg_derivatives + 3*index2+i );
  }
  //And virial
  Tensor vir( -extProduct( conn, dconn ) ); unsigned virbase = narg_derivatives + 3*getNumberOfAtoms();
  for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j ) { myvals.addDerivative( 0, virbase+3*i+j, vir(i,j) ); myvals.updateIndex( 0, virbase+3*i+j ); }
}

void TorsionsMatrix::runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return ;

  unsigned mat1s = 3*ival, ss = getPntrToArgument(1)->getShape()[1];
  unsigned nmat_ind = myvals.getNumberOfMatrixRowDerivatives();
  unsigned narg_derivatives = getPntrToArgument(0)->getNumberOfValues() + getPntrToArgument(1)->getNumberOfValues();
  std::vector<unsigned>& matrix_indices( myvals.getMatrixRowDerivativeIndices() ); unsigned ntwo_atoms = myvals.getSplitIndex();
  for(unsigned j=0; j<3; ++j) {
    matrix_indices[nmat_ind] = mat1s + j; nmat_ind++;
    matrix_indices[nmat_ind] = narg_derivatives + mat1s + j; nmat_ind++;
    for(unsigned i=1; i<ntwo_atoms; ++i) {
      unsigned ind2 = indices[i]; if( ind2>=getPntrToArgument(0)->getShape()[0] ) ind2 = indices[i] - getPntrToArgument(0)->getShape()[0];
      matrix_indices[nmat_ind] = getPntrToArgument(0)->getNumberOfStoredValues() + j*ss + ind2; nmat_ind++;
      matrix_indices[nmat_ind] = narg_derivatives + 3*indices[i] + j; nmat_ind++;
    }
  }
  unsigned base = narg_derivatives + 3*getNumberOfAtoms(); for(unsigned j=0; j<9; ++j) { matrix_indices[nmat_ind] = base + j; nmat_ind++; }
  myvals.setNumberOfMatrixRowDerivatives( nmat_ind );
}

}
}
