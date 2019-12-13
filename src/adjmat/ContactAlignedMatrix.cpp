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
#include "AlignedMatrixBase.h"
#include "core/ActionRegister.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX ALIGNED_MATRIX
/*
Adjacency matrix in which two molecule are adjacent if they are within a certain cutoff and if they have the same orientation.

As discussed in the section of the manual on \ref contactmatrix a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether
or not the \f$i\f$th and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  These matrices can then be further
analyzed using a number of other algorithms as is detailed in \cite tribello-clustering.

For this action the elements of the adjacency matrix are calculated using:

\f[
a_{ij} = \sigma_1( |\mathbf{r}_{ij}| ) \sigma_2( \mathbf{v}_i . \mathbf{v}_j )
\f]

This form of adjacency matrix can only be used if the input species are objects that lie at a point in space and that have an orientation, \f$\mathbf{v}\f$.
These orientations might represent the
orientation of a molecule, which could be calculated using \ref MOLECULES or \ref PLANES, or it might be the complex vectors calculated using the
Steinhardt parameters \ref Q3, \ref Q4 or \ref Q6.  In the expression above \f$\mathbf{r}_{ij}\f$ is the vector connecting the points in space
where objects \f$i\f$ and \f$j\f$ find themselves and \f$\sigma_1\f$ is a \ref switchingfunction that acts upon the magnitude of this vector.
\f$\sigma_2\f$ is a second \ref switchingfunction that acts on the dot product of the directors of the vectors that define the orientations of
objects \f$i\f$ and \f$j\f$.

\par Examples

The example input below is necessarily but gives you an idea of what can be achieved using this action.
The orientations and positions of four molecules are defined using the \ref MOLECULES action as the position of the
centers of mass of the two atoms specified and the direction of the vector connecting the two atoms that were specified.
A \f$4 \times 4\f$ matrix is then computed using the formula above.  The \f$ij\f$-element of this matrix tells us whether
or not atoms \f$i\f$ and \f$j\f$ are within 0.1 nm of each other and whether or not the dot-product of their orientation vectors
is greater than 0.5.  The sum of the rows of this matrix are then computed.  The sums of the \f$i\f$th row of this matrix tells us how
many of the molecules that are within the first coordination sphere of molecule \f$i\f$ have an orientation that is similar to that of
molecule \f$i\f$.  We thus calculate the number of these "coordination numbers" that are greater than 1.0 and output this quantity to a file.

\plumedfile
m1: MOLECULES MOL1=1,2 MOL2=3,4 MOL3=5,6 MOL4=7,8
mat: ALIGNED_MATRIX ATOMS=m1 SWITCH={RATIONAL R_0=0.1} ORIENTATION_SWITCH={RATIONAL R_0=0.1 D_MAX=0.5}
row: ROWSUMS MATRIX=mat MORE_THAN={RATIONAL D_0=1.0 R_0=0.1}
PRINT ARG=row.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class ContactAlignedMatrix : public AlignedMatrixBase {
private:
  Matrix<SwitchingFunction> sf;
public:
  ///
  static void registerKeywords( Keywords& keys );
  ///
  explicit ContactAlignedMatrix(const ActionOptions&);
  void readOrientationConnector( const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) override;
  double computeVectorFunction( const unsigned& iv, const unsigned& jv,
                                const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
                                Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const override;
};

PLUMED_REGISTER_ACTION(ContactAlignedMatrix,"ALIGNED_MATRIX")

void ContactAlignedMatrix::registerKeywords( Keywords& keys ) {
  AlignedMatrixBase::registerKeywords( keys );
  keys.add("numbered","ORIENTATION_SWITCH","A switching function that transforms the dot product of the input vectors.");
}

ContactAlignedMatrix::ContactAlignedMatrix( const ActionOptions& ao ):
  Action(ao),
  AlignedMatrixBase(ao)
{
  unsigned nrows, ncols, ig; retrieveTypeDimensions( nrows, ncols, ig );
  sf.resize( nrows, ncols );
  parseConnectionDescriptions("ORIENTATION_SWITCH",false,0);
}

void ContactAlignedMatrix::readOrientationConnector( const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) {
  plumed_assert( desc.size()==1 ); std::string errors; sf(j,i).set(desc[0],errors);
  if( j!=i ) sf(i,j).set(desc[0],errors);
  log.printf("  vectors in %u th and %u th groups must have a dot product that is greater than %s \n",i+1,j+1,(sf(i,j).description()).c_str() );
}

double ContactAlignedMatrix::computeVectorFunction( const unsigned& iv, const unsigned& jv,
    const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
    Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const {
  double dot_df, dot=0; dconn.zero();
  for(unsigned k=2; k<vec1.size(); ++k) dot+=vec1[k]*vec2[k];
  double f_dot = sf(iv,jv).calculate( dot, dot_df );
  for(unsigned k=2; k<vec1.size(); ++k) { dvec1[k]=dot_df*vec2[k]; dvec2[k]=dot_df*vec1[k]; }
  return f_dot;
}

}
}



