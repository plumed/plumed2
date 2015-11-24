/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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


namespace PLMD {
namespace adjmat {

class TopologyMatrix : public AdjacencyMatrixBase {
private:
/// The width to use for the kernel density estimation and the 
/// sizes of the bins to be used in kernel density estimation
  double sigma, binw, radmax2;
  std::string kerneltype;
/// The maximum number of bins that will be used 
/// This is calculated based on the dmax of the switching functions
  unsigned maxbins;
/// switching function
  Matrix<SwitchingFunction> switchingFunction;
  Matrix<SwitchingFunction> cylinder_sw;
  SwitchingFunction threshold_switch;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit TopologyMatrix(const ActionOptions&);
/// Get the number of quantities that we must compute
  unsigned getNumberOfQuantities();
/// Create the ith, ith switching function
  void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::string& desc );
/// Get the position that we should use as the center of our link cells
  Vector getLinkCellPosition( const std::vector<unsigned>& atoms ) const ;
/// This actually calculates the value of the contact function
  void calculateWeight( const unsigned& taskCode, multicolvar::AtomValuePack& myatoms ) const ;
/// This does nothing
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ;
/// Calculate the contribution from one of the atoms in the third element of the pack
  void calculateForThreeAtoms( const unsigned& iat, const Vector& d1, const double& d1_len, const double& sw,
                               const double& dfunc1, HistogramBead& bead, multicolvar::AtomValuePack& myatoms ) const ; 
///
  double transformStoredValues( const std::vector<double>& myvals, unsigned& vout, double& df ) const ;
///
  /// Used to check for connections between atoms
  bool checkForConnection( const std::vector<double>& myvals ) const ;
};

PLUMED_REGISTER_ACTION(TopologyMatrix,"TOPOLOGY_MATRIX")

void TopologyMatrix::registerKeywords( Keywords& keys ){
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("atoms","NODES","The list of atoms for which you would like to calculate the contact matrix.  The atoms involved must be specified "
                           "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
                           "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
                           "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
                           "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms","ATOMS","");
  keys.add("numbered","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.add("numbered","RADIUS","");
  keys.add("compulsory","DENSITY_THRESHOLD","");
  keys.add("compulsory","BIN_SIZE",""); keys.use("SUM");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
}

TopologyMatrix::TopologyMatrix( const ActionOptions& ao ):
Action(ao),
AdjacencyMatrixBase(ao)
{
  // Read in stuff for grid
  parse("BIN_SIZE",binw); parse("SIGMA",sigma); parse("KERNEL",kerneltype);
  // Read in threshold for density cutoff
  std::string errors, thresh_sw_str; parse("DENSITY_THRESHOLD",thresh_sw_str);
  threshold_switch.set(thresh_sw_str, errors );
  if( errors.length()>0 ) error("errors in DENSITY_THRESHOLD switching function : " + errors ); 
  log.printf("  threshold on density of atoms in cylinder equals %s\n",threshold_switch.description().c_str() );

  // The weight now does indeed have derivatives
  weightHasDerivatives=true;
  // Read in the atomic positions
  std::vector<unsigned> dims(3);
  std::vector<AtomNumber> all_atoms, atoms; parseAtomList("NODES",-1,atoms); 
  if( atoms.size()>0 ){ 
      dims[0]=dims[1]=atoms.size(); 
      for(unsigned i=0;i<atoms.size();++i) all_atoms.push_back( atoms[i] );
  } else{ dims[0]=dims[1]=colvar_label.size(); }

  switchingFunction.resize( getNumberOfNodeTypes(), getNumberOfNodeTypes() );
  parseConnectionDescriptions("SWITCH",0);
  cylinder_sw.resize( getNumberOfNodeTypes(), getNumberOfNodeTypes() );
  parseConnectionDescriptions("RADIUS",0);

  // Read in atoms 
  atoms.resize(0); parseAtomList("ATOMS",-1,atoms);
  if( atoms.size()==0 ) error("no atoms were specified");
  for(unsigned i=0;i<atoms.size();++i) all_atoms.push_back( atoms[i] );
  dims[2]=atoms.size();

  // Find the largest sf cutoff
  double sfmax=switchingFunction(0,0).get_dmax();
  double rfmax=cylinder_sw(0,0).get_dmax();
  for(unsigned i=0;i<getNumberOfNodeTypes();++i){
      for(unsigned j=0;j<getNumberOfNodeTypes();++j){
          double tsf=switchingFunction(i,j).get_dmax();
          if( tsf>sfmax ) sfmax=tsf;
          double rsf=cylinder_sw(i,j).get_dmax();
          if( rsf>rfmax ) rfmax=rsf;
      }
  }
  // Get the width of the bead
  HistogramBead bead; bead.isNotPeriodic(); 
  bead.setKernelType( kerneltype ); bead.set( 0.0, 1.0, sigma );
  double radmax = sfmax/2.0 + bead.getCutoff(); radmax2=radmax*radmax;

  // Set the link cell cutoff
  log.printf("  setting cutoffs %f %f \n",sfmax, sqrt( radmax*radmax + rfmax*rfmax ) );
  setLinkCellCutoff( sfmax, sqrt( radmax*radmax + rfmax*rfmax ) );
  // Set the maximum number of bins that we will need to compute
  maxbins = std::floor( sfmax / binw );
  
  // And request the atoms involved in this colvar
  requestAtoms( all_atoms, true, false, dims );
}

unsigned TopologyMatrix::getNumberOfQuantities(){
  return maxbins+1;
}

void TopologyMatrix::setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::string& desc ){ 
  plumed_assert( id<2 );
  if( id==0 ){
     std::string errors; switchingFunction(j,i).set(desc,errors);
     if( errors.length()!=0 ) error("problem reading switching function description " + errors);
     if( j!=i) switchingFunction(i,j).set(desc,errors);
     log.printf("  %d th and %d th multicolvar groups must be within %s\n",i+1,j+1,(switchingFunction(i,j).description()).c_str() );
  } else if( id==1 ){
     std::string errors; cylinder_sw(j,i).set(desc,errors);
     if( errors.length()!=0 ) error("problem reading switching function description " + errors);
     if( j!=i) cylinder_sw(i,j).set(desc,errors);
     log.printf("  there must be not atoms within the cylinder connections atoms of multicolvar groups %d th and %d th.  This cylinder has radius %s \n",i+1,j+1,(cylinder_sw(i,j).description()).c_str() );
  }
}

Vector TopologyMatrix::getLinkCellPosition( const std::vector<unsigned>& atoms ) const {
  Vector myatom = getPositionOfAtomForLinkCells( atoms[0] );
  return myatom + 0.5*pbcDistance( getPositionOfAtomForLinkCells( atoms[1] ), myatom );
}

void TopologyMatrix::calculateWeight( const unsigned& taskCode, multicolvar::AtomValuePack& myatoms ) const {
  Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  if( distance.modulo()<switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) ).get_dmax() ){
      myatoms.setValue(0,1);
  } else {
      myatoms.setValue(0,0);
  }
}

double TopologyMatrix::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  HistogramBead bead; bead.isNotPeriodic(); bead.setKernelType( kerneltype );

  // Initialise to zero density on all bins
  for(unsigned bin=0;bin<maxbins;++bin) myatoms.setValue(bin+1,0);
  // Calculate whether or not atoms 1 and 2 are within cutoff (can use delta here as pbc are done in atom setup)
  Vector d1 = delta( myatoms.getPosition(0), myatoms.getPosition(1) ); double d1_len = d1.modulo();
  double dfuncl, sw = switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ),
                                         getBaseColvarNumber( myatoms.getIndex(1) ) ).calculate( d1_len, dfuncl );
  d1 = d1 / d1_len;

  if( myatoms.getNumberOfAtoms()>3 ){
      for(unsigned i=2;i<myatoms.getNumberOfAtoms();++i) calculateForThreeAtoms( i, d1, d1_len, sw, dfuncl, bead, myatoms );
  } else {
      plumed_dbg_assert( myatoms.getNumberOfAtoms()==3 );
      calculateForThreeAtoms( 2, d1, d1_len, sw, dfuncl, bead, myatoms );
  }
  return myatoms.getValue(1);
}

double TopologyMatrix::transformStoredValues( const std::vector<double>& myvals, unsigned& vout, double& df  ) const {
  plumed_dbg_assert( myvals.size()==maxbins ); vout=1; double max=myvals[1];
  for(unsigned i=2;i<myvals.size();++i){
      if( myvals[i]>max ){ max=myvals[i]; vout=i; }
  }
  double ff = threshold_switch.calculate( max, df ); df*=max;
  return ff;
} 


void TopologyMatrix::calculateForThreeAtoms( const unsigned& iat, const Vector& d1, const double& d1_len, const double& sw, 
                                             const double& dfuncl, HistogramBead& bead, multicolvar::AtomValuePack& myatoms ) const {
  // Calculate if there are atoms in the cylinder (can use delta here as pbc are done in atom setup)
  Vector d2 = delta( myatoms.getPosition(0), myatoms.getPosition(iat) );
  if ( d2.modulo2()>radmax2 ) return; 
  // Now calculate projection of d2 on d1
  double proj=dotProduct(d2,d1);
  // Return if the projection is outside the length of interest
  if( proj<-bead.getCutoff() || proj>(d1_len+bead.getCutoff()) ) return;

  // Calculate the projection on the perpendicular distance from the center of the tube
  double cm = d2.modulo2() - proj*proj;

  // Now calculate the density in the cylinder
  if( cm<cylinder_sw( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(1) ) ).get_dmax2() ){
      double dfuncr, val = cylinder_sw( getBaseColvarNumber( myatoms.getIndex(0) ), 
                                        getBaseColvarNumber( myatoms.getIndex(1) ) ).calculateSqr( cm, dfuncr );

      Vector dc1, dc2, dc3, dd1, dd2, dd3;
      if( !doNotCalculateDerivatives() ){
          Tensor d1_a1; double dlen_3=d1_len*d1_len*d1_len;
          // Derivative of director connecting atom1 - atom2 wrt the position of atom 1
          d1_a1(0,0) = ( -(d1[1]*d1[1]+d1[2]*d1[2])/dlen_3 );   // dx/dx
          d1_a1(0,1) = (  d1[0]*d1[1]/dlen_3 );                 // dx/dy
          d1_a1(0,2) = (  d1[0]*d1[2]/dlen_3 );                 // dx/dz 
          d1_a1(1,0) = (  d1[1]*d1[0]/dlen_3 );                 // dy/dx
          d1_a1(1,1) = ( -(d1[0]*d1[0]+d1[2]*d1[2])/dlen_3 );   // dy/dy
          d1_a1(1,2) = (  d1[1]*d1[2]/dlen_3 );
          d1_a1(2,0) = (  d1[2]*d1[0]/dlen_3 );
          d1_a1(2,1) = (  d1[2]*d1[1]/dlen_3 );
          d1_a1(2,2) = ( -(d1[1]*d1[1]+d1[0]*d1[0])/dlen_3 ); 

          // Calculate derivatives of dot product 
          dd1 = matmul(d2, d1_a1) - matmul( Tensor::identity(), d1 );
          dd2 = matmul(d2, -d1_a1);
          dd3 = matmul( Tensor::identity(), d1 );

          // Calculate derivatives of cross product
          dc1 = dfuncr*( -d2 - proj*dd1 );
          dc2 = dfuncr*( -proj*dd2 );
          dc3 = dfuncr*( d2 - proj*dd3 );
      }

      Vector g1derivf,g2derivf,lderivf; Tensor vir;
      for(unsigned bin=0;bin<maxbins;++bin){
          bead.set( bin*binw, (bin+1)*binw, sigma ); 
          if( proj<(bin*binw-bead.getCutoff()) || proj>binw*(bin+1)+bead.getCutoff() ) continue;
          double der, contr=bead.calculateWithCutoff( proj, der );
          myatoms.addValue( 1+bin, sw*contr*val );

          if( !doNotCalculateDerivatives() ){
              g1derivf=contr*sw*dc1 + sw*val*der*dd1 - contr*val*dfuncl*d1_len*d1; addAtomDerivatives( 1+bin, 0, g1derivf, myatoms );
              g2derivf=contr*sw*dc2 + sw*val*der*dd2 + contr*val*dfuncl*d1_len*d1; addAtomDerivatives( 1+bin, 1, g2derivf, myatoms );
              lderivf=contr*sw*dc3 + sw*val*der*dd3; addAtomDerivatives( 1+bin, iat, lderivf, myatoms );
              // Virial 
              vir = -Tensor( myatoms.getPosition(0), g1derivf ) - Tensor( myatoms.getPosition(1), g2derivf ) - Tensor( myatoms.getPosition(iat), lderivf );
              myatoms.addBoxDerivatives( 1+bin, vir );
          }
      }
  }
}

bool TopologyMatrix::checkForConnection( const std::vector<double>& myvals ) const { 
  double dfake; unsigned vfake; 
  return (transformStoredValues( myvals, vfake, dfake)>epsilon);
} 

}
}

