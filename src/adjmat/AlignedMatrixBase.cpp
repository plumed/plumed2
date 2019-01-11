/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "AlignedMatrixBase.h"
#include "multicolvar/AtomValuePack.h"
#include "tools/SwitchingFunction.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace adjmat {

void AlignedMatrixBase::registerKeywords( Keywords& keys ) {
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("atoms","ATOMS","The list of molecules for which you would like to calculate the contact matrix.  The molecules involved must "
           "have an orientation so your list will be a list of the labels of \\ref mcolv or \\ref multicolvarfunction "
           "as PLUMED calculates the orientations of molecules within these operations.  Please note also that the majority "
           "of \\ref mcolv and \\ref multicolvarfunction do not calculate a molecular orientation.");
  keys.add("atoms-1","ATOMSA","The list of molecules that you would like to use for the rows of the contact matrix.  The molecules involved must "
           "have an orientation so your list will be a list of the labels of \\ref mcolv or \\ref multicolvarfunction "
           "as PLUMED calculates the orientations of molecules within these operations.  Please note also that the majority "
           "of \\ref mcolv and \\ref multicolvarfunction do not calculate a molecular orientation.");
  keys.add("atoms-1","ATOMSB","The list of molecules that you would like to use for the columns of the contact matrix.  The molecules involved must "
           "have an orientation so your list will be a list of the labels of \\ref mcolv or \\ref multicolvarfunction "
           "as PLUMED calculates the orientations of molecules within these operations.  Please note also that the majority "
           "of \\ref mcolv and \\ref multicolvarfunction do not calculate a molecular orientation.");
  keys.add("numbered","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available.");
}

AlignedMatrixBase::AlignedMatrixBase( const ActionOptions& ao ):
  Action(ao),
  AdjacencyMatrixBase(ao)
{
  // Read in the atomic positions
  readMaxTwoSpeciesMatrix( "ATOMS", "ATOMSA", "ATOMSB", true );
  unsigned nrows, ncols; retrieveTypeDimensions( nrows, ncols, ncol_t );
  if( mybasemulticolvars.size()==0 ) error("cannot use atom indices as input to this variable / input was not specified");
  if( getSizeOfInputVectors()<3 ) error("base multicolvars do not calculate an orientation");
  // Read in the switching function
  switchingFunction.resize( nrows, ncols );
  parseConnectionDescriptions("SWITCH",false,ncol_t);

  // Find the largest sf cutoff
  double sfmax=switchingFunction(0,0).get_dmax();
  for(unsigned i=0; i<getNumberOfNodeTypes(); ++i) {
    for(unsigned j=0; j<getNumberOfNodeTypes(); ++j) {
      double tsf=switchingFunction(i,j).get_dmax();
      if( tsf>sfmax ) sfmax=tsf;
    }
  }
  // And set the link cell cutoff
  setLinkCellCutoff( sfmax );
}

void AlignedMatrixBase::setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) {
  plumed_assert( id<2 );
  if( id==0 ) {
    plumed_assert( desc.size()==1 ); std::string errors; switchingFunction(j,i).set(desc[0],errors);
    if( errors.length()!=0 ) error("problem reading switching function in SWITCH keywrd description " + errors);
    if( j!=i) switchingFunction(i,j).set(desc[0],errors);
    log.printf("  %u th and %u th multicolvar groups must be within %s\n",i+1,j+1,(switchingFunction(i,j).description()).c_str() );
  } else if( id==1 ) {
    readOrientationConnector( i, j, desc );
  }
}

double AlignedMatrixBase::calculateWeight( const unsigned& taskCode, const double& weight, multicolvar::AtomValuePack& myatoms ) const {
  Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  if( distance.modulo2()<switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) )-ncol_t ).get_dmax2() ) return 1.0;
  return 0.0;
}

double AlignedMatrixBase::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  unsigned ncomp=getSizeOfInputVectors(); Vector ddistance;
  std::vector<double> orient0(ncomp), orient1(ncomp), dorient0(ncomp), dorient1(ncomp);
  Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  getInputData( 0, true, myatoms, orient0 ); getInputData( 1, true, myatoms, orient1 );
  double f_dot = computeVectorFunction( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) )-ncol_t,
                                        distance, orient0, orient1, ddistance, dorient0, dorient1 );

  // Retrieve the weight of the connection
  double dfunc, sw = switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) )-ncol_t ).calculate( distance.modulo(), dfunc );

  if( !doNotCalculateDerivatives() ) {
    addAtomDerivatives( 1, 0, (-dfunc)*f_dot*distance - sw*ddistance, myatoms );
    addAtomDerivatives( 1, 1, (+dfunc)*f_dot*distance + sw*ddistance, myatoms );
    myatoms.addBoxDerivatives( 1, (-dfunc)*f_dot*Tensor(distance,distance) - sw*extProduct( distance, ddistance ) );

    // Add derivatives of orientation
    for(unsigned k=2; k<orient0.size(); ++k) { dorient0[k]*=sw; dorient1[k]*=sw; }
    mergeInputDerivatives( 1, 2, orient0.size(), 0, dorient0, getInputDerivatives( 0, true, myatoms ), myatoms );
    mergeInputDerivatives( 1, 2, orient1.size(), 1, dorient1, getInputDerivatives( 1, true, myatoms ), myatoms );
  }
  return sw*f_dot;
}

}
}

