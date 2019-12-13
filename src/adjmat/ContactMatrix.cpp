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
#include "AdjacencyMatrixBase.h"
#include "multicolvar/AtomValuePack.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX CONTACT_MATRIX
/*
Adjacency matrix in which two atoms are adjacent if they are within a certain cutoff.

As discussed in the section of the manual on \ref contactmatrix a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether
or not the \f$i\f$th and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  These matrices can then be further
analyzed using a number of other algorithms as is detailed in \cite tribello-clustering.

For this action the elements of the contact matrix are calculated using:

\f[
 a_{ij} = \sigma( |\mathbf{r}_{ij}| )
\f]

where \f$|\mathbf{r}_{ij}|\f$ is the magnitude of the vector connecting atoms \f$i\f$ and \f$j\f$ and where \f$\sigma\f$ is a \ref switchingfunction.

\par Examples

The input shown below calculates a \f$6 \times 6\f$ matrix whose elements are equal to one if atom \f$i\f$ and atom \f$j\f$ are within 0.3 nm
of each other and which is zero otherwise.  The columns in this matrix are then summed so as to give the coordination number for each atom.
The final quantity output in the colvar file is thus the average coordination number.

\plumedfile
mat: CONTACT_MATRIX ATOMS=1-6 SWITCH={EXP D_0=0.2 R_0=0.1 D_MAX=0.66}
COLUMNSUMS MATRIX=mat MEAN LABEL=csums
PRINT ARG=csums.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace adjmat {

class ContactMatrix : public AdjacencyMatrixBase {
private:
/// Number of types that are in rows
  unsigned ncol_t;
/// switching function
  Matrix<SwitchingFunction> switchingFunction;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ContactMatrix(const ActionOptions&);
/// Create the ith, ith switching function
  void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) override;
/// This actually calculates the value of the contact function
  double calculateWeight( const unsigned& taskCode, const double& weight, multicolvar::AtomValuePack& myatoms ) const override;
/// This does nothing
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const override;
};

PLUMED_REGISTER_ACTION(ContactMatrix,"CONTACT_MATRIX")

void ContactMatrix::registerKeywords( Keywords& keys ) {
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("atoms","ATOMS","The list of atoms for which you would like to calculate the contact matrix.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("numbered","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. ");
// I added these keywords so I can test the results I get for column and row sums against output from COORDINATIONNUMBERS
/// These  should never be used in production as I think they will be much slower than COORDINATIONNUMBERS
  keys.add("hidden","ATOMSA",""); keys.add("hidden","ATOMSB","");
}

ContactMatrix::ContactMatrix( const ActionOptions& ao ):
  Action(ao),
  AdjacencyMatrixBase(ao)
{
  // Read in the atoms and setup the matrix
  readMaxTwoSpeciesMatrix( "ATOMS", "ATOMSA", "ATOMSB", true );
  unsigned nrows, ncols; retrieveTypeDimensions( nrows, ncols, ncol_t );
  switchingFunction.resize( nrows, ncols );
  // Read in the switching functions
  parseConnectionDescriptions("SWITCH",false,ncol_t);

  // Find the largest sf cutoff
  double sfmax=switchingFunction(0,0).get_dmax();
  for(unsigned i=0; i<switchingFunction.nrows(); ++i) {
    for(unsigned j=0; j<switchingFunction.ncols(); ++j) {
      double tsf=switchingFunction(i,j).get_dmax();
      if( tsf>sfmax ) sfmax=tsf;
    }
  }
  // And set the link cell cutoff
  setLinkCellCutoff( sfmax );
}

void ContactMatrix::setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) {
  plumed_assert( id==0 && desc.size()==1 ); std::string errors; switchingFunction(j,i).set(desc[0],errors);
  if( errors.length()!=0 ) error("problem reading switching function description " + errors);
  if( j!=i) switchingFunction(i,j).set(desc[0],errors);
  log.printf("  %u th and %u th multicolvar groups must be within %s\n",i+1,j+1,(switchingFunction(i,j).description()).c_str() );
}

double ContactMatrix::calculateWeight( const unsigned& taskCode, const double& weight, multicolvar::AtomValuePack& myatoms ) const {
  Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  if( distance.modulo2()<switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) - ncol_t ).get_dmax2() ) return 1.0;
  return 0.0;
}

double ContactMatrix::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  double dfunc;
  double sw = switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) - ncol_t ).calculate( distance.modulo(), dfunc );

  if( !doNotCalculateDerivatives() ) {
    addAtomDerivatives( 1, 0, (-dfunc)*distance, myatoms );
    addAtomDerivatives( 1, 1, (+dfunc)*distance, myatoms );
    myatoms.addBoxDerivatives( 1, (-dfunc)*Tensor(distance,distance) );
  }
  return sw;
}

}
}

