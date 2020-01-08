/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "tools/KernelFunctions.h"
#include "tools/Torsion.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX SMAC_MATRIX
/*
Adjacency matrix in which two molecules are adjacent if they are within a certain cutoff and if the angle between them is within certain ranges.

In this case the elements of the adjacency matrix are calculated using:

\f[
A_{ij} = \sigma(r_{ij}) \sum_n K_n(\theta_{ij})
\f]

In this expression \f$r_{ij}\f$ is the distance between molecule \f$i\f$ and molecule \f$j\f$ and \f$\sigma(r_{ij}\f$ is a
\ref switchingfunction that acts on this distance.  The $K_n functions are \ref kernelfunctions that take the torsion angle, \f$\theta_{ij}\f$, between the
internal orientation vectors for molecules \f$i\f$ and \f$j\f$ as input.  These kernel functions should be set so that they are
equal to one when the relative orientation of the molecules are as they are in the solid and equal to zero otherwise.
As the above matrix element is a product of functions it is only equal to one when the centers of mass of molecules \f$i\f$ and\f$j\f$
are with a certain distance of each other and when the molecules are aligned in some desirable way.

\par Examples

In the following example an adjacency matrix is constructed in which the \f$(i,j)\f$ element is equal to one if
molecules \f$i\f$ and \f$j\f$ are within 6 angstroms of each other and if the torsional angle between the orientations
of these molecules is close to 0 or \f$\pi\f$.  The various connected components of this matrix are determined using the
\ref DFSCLUSTERING algorithm and then the size of the largest cluster of connects molecules is output to a colvar file

\plumedfile
UNITS LENGTH=A

MOLECULES ...
MOL1=1,2,1
MOL2=5,6,5
MOL3=9,10,9
MOL4=13,14,13
MOL5=17,18,17
LABEL=m1
... MOLECULES

SMAC_MATRIX ...
   ATOMS=m1 SWITCH={RATIONAL D_0=5.99 R_0=0.1 D_MAX=6.0}
   KERNEL1={TRIANGULAR CENTER=0 SIGMA=1.0} KERNEL2={TRIANGULAR CENTER=pi SIGMA=0.6}
   LABEL=smac
... SMAC_MATRIX

dfs1: DFSCLUSTERING MATRIX=smac
cc2: CLUSTER_NATOMS CLUSTERS=dfs1 CLUSTER=1
PRINT ARG=cc2 FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class SMACMatrix : public AlignedMatrixBase {
private:
  Matrix<std::vector<KernelFunctions> > kernels;
public:
  ///
  static void registerKeywords( Keywords& keys );
  ///
  explicit SMACMatrix(const ActionOptions&);
  void readOrientationConnector( const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) override;
  double computeVectorFunction( const unsigned& iv, const unsigned& jv,
                                const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
                                Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const override;
};

PLUMED_REGISTER_ACTION(SMACMatrix,"SMAC_MATRIX")

void SMACMatrix::registerKeywords( Keywords& keys ) {
  AlignedMatrixBase::registerKeywords( keys );
  keys.add("numbered","KERNEL","The various kernels that are used to determine whether or not the molecules are aligned");
}

SMACMatrix::SMACMatrix( const ActionOptions& ao ):
  Action(ao),
  AlignedMatrixBase(ao)
{
  unsigned nrows, ncols, ig; retrieveTypeDimensions( nrows, ncols, ig );
  kernels.resize( nrows, ncols ); parseConnectionDescriptions("KERNEL",true,0);
}

void SMACMatrix::readOrientationConnector( const unsigned& iv, const unsigned& jv, const std::vector<std::string>& desc ) {
  for(int i=0; i<desc.size(); i++) {
    KernelFunctions mykernel( desc[i] );
    kernels(iv,jv).push_back( mykernel );
    if( jv!=iv ) kernels(jv,iv).push_back( mykernel );
  }
  if( kernels(iv,jv).size()==0 ) error("no kernels defined");
}

double SMACMatrix::computeVectorFunction( const unsigned& iv, const unsigned& jv,
    const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
    Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const {

  unsigned nvectors = ( vec1.size() - 2 ) / 3; plumed_assert( (vec1.size()-2)%3==0 );
  std::vector<Vector> dv1(nvectors), dv2(nvectors), tdconn(nvectors); Torsion t; std::vector<Vector> v1(nvectors), v2(nvectors);
  std::vector<std::unique_ptr<Value>> pos;
  for(unsigned i=0; i<nvectors; ++i) { pos.emplace_back( new Value() ); pos[i]->setDomain( "-pi", "pi" ); }

  for(unsigned j=0; j<nvectors; ++j) {
    for(unsigned k=0; k<3; ++k) {
      v1[j][k]=vec1[2+3*j+k]; v2[j][k]=vec2[2+3*j+k];
    }
    double angle = t.compute( v1[j], conn, v2[j], dv1[j], tdconn[j], dv2[j] );
    pos[j]->set( angle );
  }

  double ans=0; std::vector<double> deriv( nvectors ), df( nvectors, 0 );

  auto pos_ptr=Tools::unique2raw(pos);

  for(unsigned i=0; i<kernels(iv,jv).size(); ++i) {
    ans += kernels(iv,jv)[i].evaluate( pos_ptr, deriv );
    for(unsigned j=0; j<nvectors; ++j) df[j] += deriv[j];
  }
  dconn.zero(); for(unsigned j=0; j<nvectors; ++j) dconn += df[j]*tdconn[j];
  for(unsigned j=0; j<nvectors; ++j) {
    for(unsigned k=0; k<3; ++k) { dvec1[2+3*j+k]=df[j]*dv1[j][k]; dvec2[2+3*j+k]=df[j]*dv2[j][k]; }
  }
  return ans;
}

}
}



