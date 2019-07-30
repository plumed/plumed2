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
#include "tools/Matrix.h"

//+PLUMEDOC MATRIX TOPOLOGY_MATRIX
/*
Adjacency matrix in which two atoms are adjacent if they are connected topologically

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class TopologyMatrix : public AdjacencyMatrixBase {
private:
/// The width to use for the kernel density estimation and the
/// sizes of the bins to be used in kernel density estimation
  double sigma;
  std::string kerneltype;
/// The maximum number of bins that will be used
/// This is calculated based on the dmax of the switching functions
  unsigned maxbins;
/// The volume of the cells
  Matrix<double> cell_volume;
/// switching function
  Matrix<SwitchingFunction> switchingFunction;
  Matrix<SwitchingFunction> cylinder_sw;
  Matrix<SwitchingFunction> low_sf;
  double beadrad, lsfmax;
  Matrix<double> binw_mat;
  SwitchingFunction threshold_switch;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit TopologyMatrix(const ActionOptions&);
/// Get the number of quantities that we must compute
  unsigned getNumberOfQuantities() const override;
/// Create the ith, ith switching function
  void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) override;
/// This actually calculates the value of the contact function
  double calculateWeight( const unsigned& taskCode, const double& weight, multicolvar::AtomValuePack& myatoms ) const override;
/// This does nothing
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const override;
/// Calculate the contribution from one of the atoms in the third element of the pack
  void calculateForThreeAtoms( const unsigned& iat, const Vector& d1, const double& d1_len,
                               HistogramBead& bead, multicolvar::AtomValuePack& myatoms ) const ;
};

PLUMED_REGISTER_ACTION(TopologyMatrix,"TOPOLOGY_MATRIX")

void TopologyMatrix::registerKeywords( Keywords& keys ) {
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("atoms","NODES","The list of atoms for which you would like to calculate the contact matrix.  The atoms involved must be specified "
           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms","ATOMS","");
  keys.add("numbered","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available.");
  keys.add("numbered","RADIUS","");
  keys.add("numbered","CYLINDER_SWITCH","a switching function on \\f$(r_{ij}\\cdot r_{ik}-1)/r_{ij}\\f$");
  keys.add("numbered","BIN_SIZE","the size to use for the bins");
  keys.add("compulsory","DENSITY_THRESHOLD","");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.add("hidden","FAKE","");
}

TopologyMatrix::TopologyMatrix( const ActionOptions& ao ):
  Action(ao),
  AdjacencyMatrixBase(ao),
  maxbins(0)
{
  readMaxThreeSpeciesMatrix("NODES", "FAKE", "FAKE", "ATOMS", true );
  unsigned nrows, ncols, ndonor_types; retrieveTypeDimensions( nrows, ncols, ndonor_types );
  switchingFunction.resize( nrows, ncols ); parseConnectionDescriptions("SWITCH",false,ndonor_types);
  cylinder_sw.resize( nrows, ncols ); parseConnectionDescriptions("RADIUS",false,ndonor_types);
  low_sf.resize( nrows, ncols ); parseConnectionDescriptions("CYLINDER_SWITCH",false,ndonor_types);
  binw_mat.resize( nrows, ncols ); cell_volume.resize( nrows, ncols );
  parseConnectionDescriptions("BIN_SIZE",false,ndonor_types);
  // Read in stuff for grid
  parse("SIGMA",sigma); parse("KERNEL",kerneltype);
  // Read in threshold for density cutoff
  std::string errors, thresh_sw_str; parse("DENSITY_THRESHOLD",thresh_sw_str);
  threshold_switch.set(thresh_sw_str, errors );
  if( errors.length()>0 ) error("errors in DENSITY_THRESHOLD switching function : " + errors );
  log.printf("  threshold on density of atoms in cylinder equals %s\n",threshold_switch.description().c_str() );

  for(unsigned i=0; i<getNumberOfNodeTypes(); ++i) {
    for(unsigned j=0; j<getNumberOfNodeTypes(); ++j) {
      double r=cylinder_sw(i,j).get_d0() + cylinder_sw(i,j).get_r0();
      cell_volume(i,j)=binw_mat(i,j)*pi*r*r;
    }
  }

  // Find the largest sf cutoff
  lsfmax=low_sf(0,0).get_dmax();
  double sfmax=switchingFunction(0,0).get_dmax();
  double rfmax=cylinder_sw(0,0).get_dmax();
  for(unsigned i=0; i<getNumberOfNodeTypes(); ++i) {
    for(unsigned j=0; j<getNumberOfNodeTypes(); ++j) {
      double tsf=switchingFunction(i,j).get_dmax();
      if( tsf>sfmax ) sfmax=tsf;
      double rsf=cylinder_sw(i,j).get_dmax();
      if( rsf>rfmax ) rfmax=rsf;
      double lsf=low_sf(i,j).get_dmax();
      if( lsf>lsfmax ) lsfmax=lsf;
    }
  }
  // Get the width of the bead
  HistogramBead bead; bead.isNotPeriodic();
  bead.setKernelType( kerneltype ); bead.set( 0.0, 1.0, sigma );
  beadrad = bead.getCutoff();

  // Set the link cell cutoff
  log.printf("  setting cutoff %f \n", sfmax );
  setLinkCellCutoff( sfmax, std::numeric_limits<double>::max() );

  double maxsize=0;
  for(unsigned i=0; i<getNumberOfNodeTypes(); ++i) {
    for(unsigned j=0; j<getNumberOfNodeTypes(); ++j) {
      if( binw_mat(i,j)>maxsize ) maxsize=binw_mat(i,j);
    }
  }
  // Set the maximum number of bins that we will need to compute
  maxbins = std::floor( sfmax / maxsize ) + 1;
  // Need to resize functions again here to ensure that vector sizes
  // are set correctly in AdjacencyMatrixVessel
  resizeFunctions();
}

unsigned TopologyMatrix::getNumberOfQuantities() const {
  return maxbins+3;
}

void TopologyMatrix::setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) {
  plumed_assert( id<4 );
  if( id==0 ) {
    std::string errors; switchingFunction(j,i).set(desc[0],errors);
    if( errors.length()!=0 ) error("problem reading switching function description " + errors);
    if( j!=i) switchingFunction(i,j).set(desc[0],errors);
    log.printf("  %u th and %u th multicolvar groups must be within %s\n",i+1,j+1,(switchingFunction(i,j).description()).c_str() );
  } else if( id==1 ) {
    std::string errors; cylinder_sw(j,i).set(desc[0],errors);
    if( errors.length()!=0 ) error("problem reading switching function description " + errors);
    if( j!=i) cylinder_sw(i,j).set(desc[0],errors);
    log.printf("  there must be not atoms within the cylinder connections atoms of multicolvar groups %u th and %u th.  This cylinder has radius %s \n",i+1,j+1,(cylinder_sw(i,j).description()).c_str() );
  } else if( id==2 ) {
    std::string errors; low_sf(j,i).set(desc[0],errors);
    if( errors.length()!=0 ) error("problem reading switching function description " + errors);
    if( j!=i ) low_sf(i,j).set(desc[0],errors);
    log.printf("  %u th and %u th multicolvar groups must be further apart than %s\n",i+1,j+1,(low_sf(j,i).description()).c_str() );
  } else if( id==3 ) {
    Tools::convert( desc[0], binw_mat(j,i) );
    if( i!=j ) binw_mat(i,j)=binw_mat(j,i);
    log.printf("  cylinder for %u th and %u th multicolvar groups is split into bins of length %f \n",i,j,binw_mat(i,j) );
  }
}

double TopologyMatrix::calculateWeight( const unsigned& taskCode, const double& weight, multicolvar::AtomValuePack& myatoms ) const {
  Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  if( distance.modulo2()<switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) ).get_dmax2() ) return 1.0;
  return 0.0;
}

double TopologyMatrix::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  HistogramBead bead; bead.isNotPeriodic(); bead.setKernelType( kerneltype );

  // Initialise to zero density on all bins
  for(unsigned bin=0; bin<maxbins; ++bin) myatoms.setValue(bin+1,0);
  // Calculate whether or not atoms 1 and 2 are within cutoff (can use delta here as pbc are done in atom setup)
  Vector d1 = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) ); double d1_len = d1.modulo();
  d1 = d1 / d1_len;  // Convert vector into director
  AtomNumber a1 = myatoms.getAbsoluteIndex( 0 );
  AtomNumber a2 = myatoms.getAbsoluteIndex( 1 );
  for(unsigned i=2; i<myatoms.getNumberOfAtoms(); ++i) {
    AtomNumber a3 = myatoms.getAbsoluteIndex( i );
    if( a3!=a1 && a3!=a2 ) calculateForThreeAtoms( i, d1, d1_len, bead, myatoms );
  }
  // std::vector<double> binvals( 1+maxbins ); for(unsigned i=1;i<maxbins;++i) binvals[i]=myatoms.getValue(i);
  // unsigned ii; double fdf;
  //std::cout<<"HELLO DENSITY "<<myatoms.getIndex(0)<<" "<<myatoms.getIndex(1)<<" "<<transformStoredValues( binvals, ii, fdf )<<std::endl;

  // Now find the element for which the density is maximal
  unsigned vout=2; double max=myatoms.getValue( 2 );
  for(unsigned i=3; i<myatoms.getUnderlyingMultiValue().getNumberOfValues()-1; ++i) {
    if( myatoms.getValue(i)>max ) { max=myatoms.getValue(i); vout=i; }
  }
  // Calculate value and derivative of switching function between atoms 1 and 2
  double dfuncl, sw = switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ),
                                         getBaseColvarNumber( myatoms.getIndex(1) ) ).calculate( d1_len, dfuncl );
  // Transform the density
  double df, tsw = threshold_switch.calculate( max, df );
  if( !doNotCalculateDerivatives() ) {
    // Factor of d1_len is required here because d1 is normalized
    d1 *= d1_len;
    addAtomDerivatives( 2+maxbins, 0, -dfuncl*d1, myatoms );
    addAtomDerivatives( 2+maxbins, 1, dfuncl*d1, myatoms );
    myatoms.addBoxDerivatives( 2+maxbins, (-dfuncl)*Tensor(d1,d1) );
    // Update active atoms so that next bit works
    updateActiveAtoms( myatoms );
    // Now finish caclulation of derivatives
    MultiValue& myvals=myatoms.getUnderlyingMultiValue();
    for(unsigned jd=0; jd<myvals.getNumberActive(); ++jd) {
      unsigned ider=myvals.getActiveIndex(jd);
      myvals.addDerivative( 1, ider, sw*df*max*myvals.getDerivative( vout, ider ) + tsw*myvals.getDerivative( 2+maxbins, ider ) );
    }
  }
  return sw*tsw;
}

void TopologyMatrix::calculateForThreeAtoms( const unsigned& iat, const Vector& d1, const double& d1_len,
    HistogramBead& bead, multicolvar::AtomValuePack& myatoms ) const {
  // Calculate if there are atoms in the cylinder (can use delta here as pbc are done in atom setup)
  Vector d2 = getSeparation( myatoms.getPosition(0), myatoms.getPosition(iat) );
  // Now calculate projection of d2 on d1
  double proj=dotProduct(d2,d1);
  // This tells us if we are outside the end of the cylinder
  double excess = proj - d1_len;
  // Return if we are outside of the cylinder as calculated based on excess
  if( excess>low_sf( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) ).get_dmax() ) return;
  // Find the length of the cylinder
  double binw = binw_mat( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) );
  double lcylinder = (std::floor( d1_len / binw ) + 1)*binw;
  // Return if the projection is outside the length of interest
  if( proj<-bead.getCutoff() || proj>(lcylinder+bead.getCutoff()) ) return;

  // Calculate the excess swiching function
  double edf, eval = low_sf( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) ).calculate( excess, edf );
  // Calculate the projection on the perpendicular distance from the center of the tube
  double cm = d2.modulo2() - proj*proj;

  // Now calculate the density in the cylinder
  if( cm<cylinder_sw( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) ).get_dmax2() ) {
    double dfuncr, val = cylinder_sw( getBaseColvarNumber( myatoms.getIndex(0) ),
                                      getBaseColvarNumber( myatoms.getIndex(1) ) ).calculateSqr( cm, dfuncr );
    double cellv = cell_volume( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) );
    Vector dc1, dc2, dc3, dd1, dd2, dd3, de1, de2, de3;
    if( !doNotCalculateDerivatives() ) {
      Tensor d1_a1;
      // Derivative of director connecting atom1 - atom2 wrt the position of atom 1
      d1_a1(0,0) = ( -(d1[1]*d1[1]+d1[2]*d1[2])/d1_len );   // dx/dx
      d1_a1(0,1) = (  d1[0]*d1[1]/d1_len );                 // dx/dy
      d1_a1(0,2) = (  d1[0]*d1[2]/d1_len );                 // dx/dz
      d1_a1(1,0) = (  d1[1]*d1[0]/d1_len );                 // dy/dx
      d1_a1(1,1) = ( -(d1[0]*d1[0]+d1[2]*d1[2])/d1_len );   // dy/dy
      d1_a1(1,2) = (  d1[1]*d1[2]/d1_len );
      d1_a1(2,0) = (  d1[2]*d1[0]/d1_len );
      d1_a1(2,1) = (  d1[2]*d1[1]/d1_len );
      d1_a1(2,2) = ( -(d1[1]*d1[1]+d1[0]*d1[0])/d1_len );

      // Calculate derivatives of dot product
      dd1 = matmul(d2, d1_a1) - d1;
      dd2 = matmul(d2, -d1_a1);
      dd3 = d1;

      // Calculate derivatives of cross product
      dc1 = dfuncr*( -d2 - proj*dd1 );
      dc2 = dfuncr*( -proj*dd2 );
      dc3 = dfuncr*( d2 - proj*dd3 );

      // Calculate derivatives of excess
      de1 = edf*excess*( dd1 + d1 );
      de2 = edf*excess*( dd2 - d1 );
      de3 = edf*excess*dd3;
    }

    Vector pos1 = myatoms.getPosition(0) + d1_len*d1;
    Vector pos2 = myatoms.getPosition(0) + d2;
    Vector g1derivf,g2derivf,lderivf; Tensor vir;
    for(unsigned bin=0; bin<maxbins; ++bin) {
      bead.set( bin*binw, (bin+1)*binw, sigma );
      if( proj<(bin*binw-bead.getCutoff()) || proj>binw*(bin+1)+bead.getCutoff() ) continue;
      double der, contr=bead.calculateWithCutoff( proj, der ) / cellv; der /= cellv;
      myatoms.addValue( 2+bin, contr*val*eval );

      if( !doNotCalculateDerivatives() ) {
        g1derivf=contr*eval*dc1 + val*eval*der*dd1 + contr*val*de1;
        addAtomDerivatives( 2+bin, 0, g1derivf, myatoms );
        g2derivf=contr*eval*dc2 + val*eval*der*dd2 + contr*val*de2;
        addAtomDerivatives( 2+bin, 1, g2derivf, myatoms );
        lderivf=contr*eval*dc3 + val*eval*der*dd3 + contr*val*de3;
        addAtomDerivatives( 2+bin, iat, lderivf, myatoms );
        // Virial
        vir = -Tensor( myatoms.getPosition(0), g1derivf ) - Tensor( pos1, g2derivf ) - Tensor( pos2, lderivf );
        myatoms.addBoxDerivatives( 2+bin, vir );
      }
    }
  }
}

}
}

