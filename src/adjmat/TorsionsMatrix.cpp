/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "VectorProductMatrix.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionRegister.h"
#include "tools/Torsion.h"

namespace PLMD {
namespace adjmat {

class TorsionsMatrix : public VectorProductMatrix {
private:
public:
  static void registerKeywords( Keywords& keys );
  explicit TorsionsMatrix(const ActionOptions&);
  double computeVectorProduct( const unsigned& index1, const unsigned& index2,
                               const std::vector<double>& vec1, const std::vector<double>& vec2,
                               std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(TorsionsMatrix,"TORSIONS_MATRIX")

void TorsionsMatrix::registerKeywords( Keywords& keys ) {
  VectorProductMatrix::registerKeywords( keys );
  keys.add("atoms","POSITIONS","the positions to use for the molecules");
  keys.add("atoms-1","POSITIONSA","the positions to use for the molecules specified using GROUPA");
  keys.add("atoms-1","POSITIONSB","the positions to use for the molecules specified using GROUPB");
}

TorsionsMatrix::TorsionsMatrix(const ActionOptions& ao):
  Action(ao),
  VectorProductMatrix(ao)
{
  setPeriodic( "-pi", "pi" );  // Make sure that the periodicity of the value is set
  unsigned nargs=getNumberOfArguments(); if( ncol_args>0 ) nargs /= 2;
  if( ncol_args>0 ) {
    std::vector<AtomNumber> atoms_a; parseAtomList("POSITIONSA", atoms_a );
    if( atoms_a.size()!=getPntrToArgument(0)->getShape()[0] ) error("mismatch between number of atoms specified using POSITIONSA and number of arguments in vector input");
    log.printf("  using positions of these atoms for vectors in rows of matrix \n");
    for(unsigned int i=0; i<atoms_a.size(); ++i) {
      if ( (i+1) % 25 == 0 ) log.printf("  \n");
      log.printf("  %d", atoms_a[i].serial());
    }
    log.printf("\n"); std::vector<AtomNumber> atoms_b; parseAtomList("POSITIONSB", atoms_b );
    if( atoms_b.size()!=getPntrToArgument(ncol_args)->getShape()[0] ) error("mismatch between number of atoms specified using POSITIONSB and number of arguments in vector input");
    log.printf("  using positions of these atoms for vectors in columns of matrix \n");
    for(unsigned i=0; i<atoms_b.size(); ++i) {
      if ( (i+1) % 25 == 0 ) log.printf("  \n");
      log.printf("  %d", atoms_b[i].serial()); atoms_a.push_back( atoms_b[i] );
    }
    log.printf("\n"); std::vector<Value*> args( getArguments() ); requestAtoms( atoms_a ); requestArguments( args, false );
  } else {
    std::vector<AtomNumber> atoms; parseAtomList("POSITIONS", atoms );
    if( atoms.size()!=getPntrToArgument(0)->getShape()[0] ) error("mismatch between number of atoms specified and number of arguments in vector");
    log.printf("  using positions of these atoms when calculating torsions \n");
    for(unsigned int i=0; i<atoms.size(); ++i) {
      if ( (i+1) % 25 == 0 ) log.printf("  \n");
      log.printf("  %d", atoms[i].serial());
    }
    log.printf("\n"); std::vector<Value*> args( getArguments() ); requestAtoms( atoms ); requestArguments( args, false );
  }
  forcesToApply.resize( getNumberOfDerivatives() );
}

double TorsionsMatrix::computeVectorProduct( const unsigned& index1, const unsigned& index2,
    const std::vector<double>& vec1, const std::vector<double>& vec2,
    std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const {
  plumed_dbg_assert( vec1.size()==3 && vec2.size()==3 );
  Vector v1, v2, dv1, dv2, dconn;
  // Compute the distance connecting the two centers
  Vector conn=pbcDistance( getPosition(index1), getPosition(index2) );
  for(unsigned i=0; i<3; ++i) { v1[i]=vec1[i]; v2[i]=vec2[i]; }
  // Evaluate angle
  Torsion t; double angle = t.compute( v1, conn, v2, dv1, dconn, dv2 );
  if( !doNotCalculateDerivatives() ) {
    unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
    for(unsigned i=0; i<3; ++i) { dvec1[i]=dv1[i]; dvec2[i]=dv2[i]; }
    myvals.addDerivative( ostrn, narg_derivatives + 3*index1+0, -dconn[0] ); myvals.addDerivative( ostrn, narg_derivatives + 3*index2+0, dconn[0] );
    myvals.updateIndex( ostrn, narg_derivatives + 3*index1+0 ); myvals.updateIndex( ostrn, narg_derivatives + 3*index2+0 );
    myvals.addDerivative( ostrn, narg_derivatives + 3*index1+1, -dconn[1] ); myvals.addDerivative( ostrn, narg_derivatives + 3*index2+1, dconn[1] );
    myvals.updateIndex( ostrn, narg_derivatives + 3*index1+1 ); myvals.updateIndex( ostrn, narg_derivatives + 3*index2+1 );
    myvals.addDerivative( ostrn, narg_derivatives + 3*index1+2, -dconn[2] ); myvals.addDerivative( ostrn, narg_derivatives + 3*index2+2, dconn[2] );
    myvals.updateIndex( ostrn, narg_derivatives + 3*index1+2 ); myvals.updateIndex( ostrn, narg_derivatives + 3*index2+2 );
    //And virial
    Tensor vir( -extProduct( conn, dconn ) ); unsigned virbase = narg_derivatives + 3*getNumberOfAtoms();
    for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j ) {
        myvals.addDerivative( ostrn, virbase+3*i+j, vir(i,j) ); myvals.updateIndex( ostrn, virbase+3*i+j );
      }
  }
  return angle;
}

}
}
