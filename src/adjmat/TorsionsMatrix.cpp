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
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit TorsionsMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives() const ;
  double computeVectorProduct( const unsigned& index1, const unsigned& index2,
                               const std::vector<double>& vec1, const std::vector<double>& vec2,
                               std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(TorsionsMatrix,"TORSIONS_MATRIX")
PLUMED_REGISTER_SHORTCUT(TorsionsMatrix,"SMAC")

void TorsionsMatrix::shortcutKeywords( Keywords& keys ) {
  keys.add("optional","SPECIES","");
  keys.add("optional","SPECIESA","");
  keys.add("optional","SPECIESB","");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.add("numbered","KERNEL","The kernels used in the function of the angle");
  keys.add("optional","SWITCH_COORD","This keyword is used to define the coordination switching function.");
  keys.reset_style("KERNEL","optional");
  multicolvar::MultiColvarBase::shortcutKeywords( keys );
}

void TorsionsMatrix::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                     const std::map<std::string,std::string>& keys,
                                     std::vector<std::vector<std::string> >& actions ) {
  // Create the matrices
  std::vector<std::string> cmap_input; cmap_input.push_back(lab + "_cmap:"); cmap_input.push_back("CONTACT_MATRIX");
  std::vector<std::string> tpmat_input; tpmat_input.push_back(lab + "_tpmat:"); tpmat_input.push_back("TORSIONS_MATRIX");
  if( keys.count("SPECIES") ) {
    std::string sp_lab = keys.find("SPECIES")->second;
    tpmat_input.push_back("GROUP1=" + sp_lab + ".x");
    tpmat_input.push_back("GROUP2=" + sp_lab + ".y");
    tpmat_input.push_back("GROUP3=" + sp_lab + ".z");
    tpmat_input.push_back("POSITIONS=" + sp_lab ); cmap_input.push_back("GROUP=" + sp_lab );
  } else if( keys.count("SPECIESA") ) {
    std::string sp_laba = keys.find("SPECIESA")->second; std::string sp_labb = keys.find("SPECIESB")->second;
    cmap_input.push_back( "GROUPA=" + sp_laba ); cmap_input.push_back( "GROUPB=" + sp_labb );
    tpmat_input.push_back( "GROUPA1=" + sp_laba + ".x"); tpmat_input.push_back( "GROUPB1=" + sp_labb + ".x");
    tpmat_input.push_back( "GROUPA2=" + sp_laba + ".y"); tpmat_input.push_back( "GROUPB2=" + sp_labb + ".y");
    tpmat_input.push_back( "GROUPA3=" + sp_laba + ".z"); tpmat_input.push_back( "GROUPB3=" + sp_labb + ".z");
    tpmat_input.push_back("POSITIONSA=" + sp_laba ); tpmat_input.push_back("POSITIONSB=" + sp_labb );
  }
  cmap_input.push_back( "SWITCH=" + keys.find("SWITCH")->second );
  actions.push_back( cmap_input ); actions.push_back( tpmat_input );
  // Now need the Gaussians
  std::vector<std::string> kmap_input; kmap_input.push_back(lab + "_ksum:");
  kmap_input.push_back("COMBINE"); kmap_input.push_back("PERIODIC=NO");
  for(unsigned i=1;; ++i) {
    std::string istr; Tools::convert( i, istr );
    if( !keys.count("KERNEL" + istr ) ) { break; }
    std::vector<std::string> input; input.push_back( lab + "_kf" + istr + ":" ); input.push_back("KERNEL");
    input.push_back("ARG1=" + lab + "_tpmat"); input.push_back("KERNEL=" + keys.find("KERNEL" + istr)->second );
    actions.push_back( input ); kmap_input.push_back("ARG" + istr + "=" + lab + "_kf" + istr );
  }
  actions.push_back( kmap_input );
  // Now create the product matrix
  std::vector<std::string> pmat_input; pmat_input.push_back(lab + "_prod:");
  pmat_input.push_back("MATHEVAL"); pmat_input.push_back("ARG1=" + lab + "_cmap.w");
  pmat_input.push_back("ARG2=" + lab + "_ksum"); pmat_input.push_back("FUNC=x*y");
  pmat_input.push_back("PERIODIC=NO"); actions.push_back( pmat_input );
  // Now the sum of coordination numbers times the switching functions
  std::vector<std::string> coord_input_numer; coord_input_numer.push_back(lab + ":");
  coord_input_numer.push_back("COORDINATIONNUMBER"); coord_input_numer.push_back("WEIGHT=" + lab + "_prod");
  actions.push_back( coord_input_numer );
  // And just the sum of the coordination numbers
  std::vector<std::string> coord_input_denom; coord_input_denom.push_back(lab + "_denom:");
  coord_input_denom.push_back("COORDINATIONNUMBER"); coord_input_denom.push_back("WEIGHT=" + lab + "_cmap.w");
  actions.push_back( coord_input_denom );
  // And the transformed switching functions
  std::vector<std::string> coord_lessthan; coord_lessthan.push_back(lab + "_mtdenom:");
  coord_lessthan.push_back("MORE_THAN"); coord_lessthan.push_back("ARG1=" + lab + "_denom");
  coord_lessthan.push_back("SWITCH=" + keys.find("SWITCH_COORD")->second );
  actions.push_back( coord_lessthan );
  // And matheval to get the final quantity
  std::vector<std::string> matheval_input; matheval_input.push_back(lab + "_smac:");
  matheval_input.push_back("MATHEVAL");
  matheval_input.push_back("ARG1=" + lab );
  matheval_input.push_back("ARG2=" + lab + "_mtdenom");
  matheval_input.push_back("ARG3=" + lab + "_denom");
  matheval_input.push_back("FUNC=(x*y)/z");
  matheval_input.push_back("PERIODIC=NO"); actions.push_back( matheval_input );
  // And this expands everything
  multicolvar::MultiColvarBase::expandFunctions( lab, lab + "_smac", "", words, keys, actions );
}

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

unsigned TorsionsMatrix::getNumberOfDerivatives() const  {
  unsigned nat_der = 3*getNumberOfAtoms()+9; plumed_assert( getPntrToArgument(0)->getRank()!=0 );
  if( ncol_args>0 ) return nat_der + (getPntrToArgument(0)->getShape()[0]+getPntrToArgument(ncol_args)->getShape()[0])*getNumberOfArguments()/2;
  return nat_der + getPntrToArgument(0)->getShape()[0]*getNumberOfArguments();
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
