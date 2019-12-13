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
#include "tools/HistogramBead.h"
#include "tools/Angle.h"
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX HBOND_MATRIX
/*
Adjacency matrix in which two atoms are adjacent if there is a hydrogen bond between them.

As discussed in the section of the manual on \ref contactmatrix a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether
or not the \f$i\f$th and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  These matrices can then be further
analyzed using a number of other algorithms as is detailed in \cite tribello-clustering.

For this action the elements of the adjacency matrix are calculated using:

\f[
a_{ij} = \sigma_{oo}( |\mathbf{r}_{ij}| ) \sum_{k=1}^N \sigma_{oh}( |\mathbf{r}_{ik}| ) \sigma_{\theta}( \theta_{kij} )
\f]

This expression was derived by thinking about how to detect if there is a hydrogen bond between atoms \f$i\f$ and \f$j\f$.  The notion is that
if the hydrogen bond is present atoms \f$i\f$ and \f$j\f$ should be within a certain cutoff distance.  In addition, there should be a hydrogen
within a certain cutoff distance of atom \f$i\f$ and this hydrogen should lie on or close to the vector connecting atoms \f$i\f$ and \f$j\f$.
As such \f$\sigma_{oo}( |\mathbf{r}_{ij}| )\f$ is a \ref switchingfunction that acts on the modulus of the vector connecting atom \f$i\f$ to atom
\f$j\f$.  The sum over \f$k\f$ then runs over all the hydrogen atoms that are specified using using HYDROGEN keyword.  \f$\sigma_{oh}(|\mathbf{r}_{ik}|)\f$
is a \ref switchingfunction that acts on the modulus of the vector connecting atom \f$i\f$ to atom \f$k\f$ and \f$\sigma_{\theta}(\theta_{kij})\f$
is a \ref switchingfunction that acts on the angle between the vector connecting atoms \f$i\f$ and \f$j\f$ and the vector connecting atoms \f$i\f$ and
\f$k\f$.

It is important to note that hydrogen bonds, unlike regular bonds, are asymmetric. In other words, the hydrogen atom does not sit at the
mid point between the two other atoms in this three-center bond.  As a result of this adjacency matrices calculated using \ref HBOND_MATRIX are not
symmetric like those calculated by \ref CONTACT_MATRIX.  One consequence of this fact is that the quantities found by performing \ref ROWSUMS and
\ref COLUMNSUMS on a square \ref HBOND_MATRIX are not the same as they would be if you performed \ref ROWSUMS and
\ref COLUMNSUMS on a square \ref CONTACT_MATRIX.

\par Examples

The following input can be used to analyze the number of hydrogen bonds each of the oxygen atoms in a box of water participates in.  Each
water molecule can participate in a hydrogen bond in one of two ways.  It can either donate one of its hydrogen atom to the neighboring oxygen or
it can accept a bond between the hydrogen of a neighboring water molecule and its own oxygen.  The input below allows you to output information
on the number of hydrogen bonds each of the water molecules donates and accepts.  This information is output in two xyz files which each contain
five columns of data.  The first four of these columns are a label for the atom and the x, y and z position of the oxygen.  The last column is then
the number of accepted/donated hydrogen bonds.

\plumedfile
mat: HBOND_MATRIX ATOMS=1-192:3 HYDROGENS=2-192:3,3-192:3 SWITCH={RATIONAL R_0=3.20} HSWITCH={RATIONAL R_0=2.30} ASWITCH={RATIONAL R_0=0.167pi} SUM
rsums: ROWSUMS MATRIX=mat MEAN
csums: COLUMNSUMS MATRIX=mat MEAN
DUMPMULTICOLVAR DATA=rsums FILE=donors.xyz
DUMPMULTICOLVAR DATA=csums FILE=acceptors.xyz
\endplumedfile

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace adjmat {

class HBondMatrix : public AdjacencyMatrixBase {
private:
  unsigned ndonor_types;
/// switching function
  Matrix<SwitchingFunction> distanceOOSwitch;
  Matrix<SwitchingFunction> distanceOHSwitch;
  Matrix<SwitchingFunction> angleSwitch;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit HBondMatrix(const ActionOptions&);
/// Create the ith, ith switching function
  void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) override;
/// This actually calculates the value of the contact function
  double calculateWeight( const unsigned& taskCode, const double& weight, multicolvar::AtomValuePack& myatoms ) const override;
/// This does nothing
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const override;
///
  double calculateForThree( const unsigned& iat, const unsigned& ano, const unsigned& dno, const Vector& ood,
                            const double& ood_df, const double& ood_sw, multicolvar::AtomValuePack& myatoms ) const;
};

PLUMED_REGISTER_ACTION(HBondMatrix,"HBOND_MATRIX")

void HBondMatrix::registerKeywords( Keywords& keys ) {
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("atoms","ATOMS","The list of atoms which can be part of a hydrogen bond.  When this command is used the set of atoms that can donate a "
           "hydrogen bond is assumed to be the same as the set of atoms that can form hydrogen bonds.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms","HYDROGENS","The list of hydrogen atoms that can form part of a hydrogen bond.  The atoms must be specified using a comma separated list, "
           "an index range or by using a \\ref GROUP.  A list of hydrogen atoms is always required even if you specify the other atoms using "
           "DONORS and ACCEPTORS as described below.");
  keys.add("atoms-2","DONORS","The list of atoms which can donate a hydrogen bond.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms-2","ACCEPTORS","The list of atoms which can accept a hydrogen bond.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("numbered","SWITCH","The \\ref switchingfunction that specifies how close a pair of atoms must be together for there to be a hydrogen bond between them");
  keys.add("numbered","HSWITCH","The \\ref switchingfunction that specifies how close the hydrogen must be to the donor atom of the hydrogen bond for it to be "
           "considered a hydrogen bond");
  keys.add("numbered","ASWITCH","A \\ref switchingfunction that is used to specify what the angle between the vector connecting the donor atom to the acceptor atom and "
           "the vector connecting the donor atom to the hydrogen must be in order for it considered to be a hydrogen bond");
  keys.use("SUM");
}

HBondMatrix::HBondMatrix( const ActionOptions& ao ):
  Action(ao),
  AdjacencyMatrixBase(ao)
{
  readMaxThreeSpeciesMatrix( "ATOMS", "DONORS", "ACCEPTORS", "HYDROGENS", false );
  unsigned nrows, ncols; retrieveTypeDimensions( nrows, ncols, ndonor_types );
  distanceOOSwitch.resize( nrows, ncols ); distanceOHSwitch.resize( nrows, ncols ); angleSwitch.resize( nrows, ncols );
  parseConnectionDescriptions("SWITCH",false,ndonor_types);
  parseConnectionDescriptions("HSWITCH",false,ndonor_types);
  parseConnectionDescriptions("ASWITCH",false,ndonor_types);

  // Find the largest sf cutoff
  double sfmax=distanceOOSwitch(0,0).get_dmax();
  for(unsigned i=0; i<getNumberOfNodeTypes(); ++i) {
    for(unsigned j=0; j<getNumberOfNodeTypes(); ++j) {
      double tsf=distanceOOSwitch(i,j).get_dmax();
      if( tsf>sfmax ) sfmax=tsf;
    }
  }
  // Set the link cell cutoff
  setLinkCellCutoff( sfmax );
}

void HBondMatrix::setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) {
  plumed_assert( id<3 && desc.size()==1 );
  if( id==0 ) {
    std::string errors; distanceOOSwitch(j,i).set(desc[0],errors);
    if( errors.length()!=0 ) error("problem reading switching function description " + errors);
    if( j!=i) distanceOOSwitch(i,j).set(desc[0],errors);
    log.printf("  atoms of type %u and %u must be within %s\n",i+1,j+1,(distanceOOSwitch(i,j).description()).c_str() );
  } else if( id==1 ) {
    std::string errors; distanceOHSwitch(j,i).set(desc[0],errors);
    if( errors.length()!=0 ) error("problem reading switching function description " + errors);
    if( j!=i) distanceOHSwitch(i,j).set(desc[0],errors);
    log.printf("  for atoms of type %u and %u the OH distance must be less than %s \n",i+1,j+1,(distanceOHSwitch(i,j).description()).c_str() );
  } else if( id==2 ) {
    std::string errors; angleSwitch(j,i).set(desc[0],errors);
    if( errors.length()!=0 ) error("problem reading switching function description " + errors);
    if( j!=i) angleSwitch(i,j).set(desc[0],errors);
    log.printf("  for atoms of type %u and %u the OOH angle must be less than %s \n",i+1,j+1,(angleSwitch(i,j).description()).c_str() );
  }
}

double HBondMatrix::calculateWeight( const unsigned& taskCode, const double& weight, multicolvar::AtomValuePack& myatoms ) const {
  // Ensure we skip diagonal elements of square matrix
  if( myatoms.getIndex(0)==myatoms.getIndex(1) ) return 0.0;

  Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  if( distance.modulo2()<distanceOOSwitch( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) ).get_dmax2() ) return 1.0;
  return 0.0;
}

double HBondMatrix::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  Vector ood = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) ); double ood_l = ood.modulo(); // acceptor - donor
  double ood_df, ood_sw=distanceOOSwitch( getBaseColvarNumber( myatoms.getIndex(0) ),
                                          getBaseColvarNumber( myatoms.getIndex(1) ) ).calculate( ood_l, ood_df );

  // Get the base colvar numbers
  unsigned ano, dno = getBaseColvarNumber( myatoms.getIndex(0) );
  if( ndonor_types==0 ) ano = getBaseColvarNumber( myatoms.getIndex(1) );
  else ano = getBaseColvarNumber( myatoms.getIndex(1) ) - ndonor_types;

  double value=0;
  if( myatoms.getNumberOfAtoms()>3 ) {
    for(unsigned i=2; i<myatoms.getNumberOfAtoms(); ++i) value+=calculateForThree( i, ano, dno, ood, ood_df, ood_sw,  myatoms );
  } else {
    plumed_dbg_assert( myatoms.getNumberOfAtoms()==3 );
    value=calculateForThree( 2, ano, dno, ood, ood_df, ood_sw, myatoms );
  }
  return value;
}

double HBondMatrix::calculateForThree( const unsigned& iat, const unsigned& ano, const unsigned& dno, const Vector& ood,
                                       const double& ood_df, const double& ood_sw, multicolvar::AtomValuePack& myatoms ) const {
  Vector ohd=getSeparation( myatoms.getPosition(0), myatoms.getPosition(iat) ); double ohd_l=ohd.modulo();
  double ohd_df, ohd_sw=distanceOHSwitch( getBaseColvarNumber( myatoms.getIndex(0) ),
                                          getBaseColvarNumber( myatoms.getIndex(1) ) ).calculate( ohd_l, ohd_df );

  Angle a; Vector ood_adf, ohd_adf; double angle=a.compute( ood, ohd, ood_adf, ohd_adf );
  double angle_df, angle_sw=angleSwitch( getBaseColvarNumber( myatoms.getIndex(0) ),
                                         getBaseColvarNumber( myatoms.getIndex(1) ) ).calculate( angle, angle_df );

  if( !doNotCalculateDerivatives() ) {
    addAtomDerivatives( 1, 0, angle_sw*ohd_sw*(-ood_df)*ood + angle_sw*ood_sw*(-ohd_df)*ohd + ood_sw*ohd_sw*angle_df*angle*(-ood_adf-ohd_adf), myatoms );
    addAtomDerivatives( 1, 1, angle_sw*ohd_sw*(+ood_df)*ood + ood_sw*ohd_sw*angle_df*angle*ood_adf, myatoms );
    addAtomDerivatives( 1, iat, angle_sw*ood_sw*(+ohd_df)*ohd + ood_sw*ohd_sw*angle_df*angle*ohd_adf, myatoms );
    myatoms.addBoxDerivatives( 1, angle_sw*ohd_sw*(-ood_df)*Tensor(ood,ood) + angle_sw*ood_sw*(-ohd_df)*Tensor(ohd,ohd)
                               -ood_sw*ohd_sw*angle_df*angle*(Tensor(ood,ood_adf)+Tensor(ohd,ohd_adf)) );
  }
  return ood_sw*ohd_sw*angle_sw;
}

}
}

