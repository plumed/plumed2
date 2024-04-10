/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "tools/SwitchingFunction.h"
#include "tools/HistogramBead.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

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
  double cell_volume;
/// switching function
  SwitchingFunction switchingFunction;
  SwitchingFunction cylinder_sw;
  SwitchingFunction low_sf;
  double binw_mat;
  SwitchingFunction threshold_switch;
public:
  static void registerKeywords( Keywords& keys );
  explicit TopologyMatrix(const ActionOptions&);
// active methods:
  double calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const override;
};

PLUMED_REGISTER_ACTION(TopologyMatrix,"TOPOLOGY_MATRIX")

void TopologyMatrix::registerKeywords( Keywords& keys ) {
  AdjacencyMatrixBase::registerKeywords( keys );
  keys.add("atoms","BACKGROUND_ATOMS","the list of atoms that should be considered as part of the background density");
  keys.add("compulsory","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.add("compulsory","RADIUS","");
  keys.add("compulsory","CYLINDER_SWITCH","a switching function on ( r_ij . r_ik - 1 )/r_ij");
  keys.add("compulsory","BIN_SIZE","the size to use for the bins");
  keys.add("compulsory","DENSITY_THRESHOLD","");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
}

TopologyMatrix::TopologyMatrix(const ActionOptions&ao):
  Action(ao),
  AdjacencyMatrixBase(ao)
{
  std::string sfinput,errors; parse("SWITCH",sfinput);
  if( sfinput.length()==0 ) error("could not find SWITCH keyword");
  switchingFunction.set(sfinput,errors);
  if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );

  std::string hsfinput; parse("CYLINDER_SWITCH",hsfinput);
  if( hsfinput.length()==0 ) error("could not find CYLINDER_SWITCH keyword");
  low_sf.set(hsfinput,errors);
  if( errors.length()!=0 ) error("problem reading CYLINDER_SWITCH keyword : " + errors );

  std::string asfinput; parse("RADIUS",asfinput);
  if( asfinput.length()==0 ) error("could not find RADIUS keyword");
  cylinder_sw.set(asfinput,errors);
  if( errors.length()!=0 ) error("problem reading RADIUS keyword : " + errors );

  std::string tsfinput; parse("DENSITY_THRESHOLD",tsfinput);
  if( tsfinput.length()==0 ) error("could not find DENSITY_THRESHOLD keyword");
  threshold_switch.set(tsfinput,errors);
  if( errors.length()!=0 ) error("problem reading DENSITY_THRESHOLD keyword : " + errors );
  // Read in stuff for grid
  parse("SIGMA",sigma); parse("KERNEL",kerneltype); parse("BIN_SIZE",binw_mat);

  // Set the link cell cutoff
  setLinkCellCutoff( true, switchingFunction.get_dmax(), std::numeric_limits<double>::max() );
  // Set the number of bins
  maxbins = std::floor( switchingFunction.get_dmax() / binw_mat ) + 1;
  // Set the cell volume
  double r=cylinder_sw.get_d0() + cylinder_sw.get_r0();
  cell_volume=binw_mat*pi*r*r;

  // And check everything has been read in correctly
  checkRead();
}

double TopologyMatrix::calculateWeight( const Vector& pos1, const Vector& pos2, const unsigned& natoms, MultiValue& myvals ) const {
  // Compute switching function on distance between atoms
  Vector distance = pbcDistance( pos1, pos2 ); double len2 = distance.modulo2();
  if( len2>switchingFunction.get_dmax2() ) return 0.0;
  double dfuncl, sw = switchingFunction.calculateSqr( len2, dfuncl );

  // Now run through all sea atoms
  HistogramBead bead; bead.isNotPeriodic(); bead.setKernelType( kerneltype );
  Vector g1derivf,g2derivf,lderivf; Tensor vir; double binlength = maxbins * binw_mat;
  MultiValue tvals( maxbins, myvals.getNumberOfDerivatives() );
  for(unsigned i=0; i<natoms; ++i) {
    // Position of sea atom (this will be the origin)
    Vector d2 = getPosition(i,myvals);
    // Vector connecting sea atom and first in bond taking pbc into account
    Vector d20 = pbcDistance( d2, pos1 );
    // Vector connecting sea atom and second in bond taking pbc into account
    Vector d21 = pbcDistance( d2, pos2 );
    // Now length of bond modulus and so on -- no pbc here as we want sea atom in middle
    Vector d1 = delta( d20, d21 ); double d1_len = d1.modulo(); d1 = d1 / d1_len;
    // Switching function on distance between nodes
    if( d1_len>switchingFunction.get_dmax() ) continue ;
    // Ensure that the center of the bins are on the center of the bond connecting the two atoms
    double start2atom = 0.5*(binlength-d1_len); Vector dstart = d20 - start2atom*d1;
    // Now calculate projection of axis of cylinder
    double proj=dotProduct(-dstart,d1);
    // Calculate length of vector connecting start of cylinder to first atom
    // And now work out projection on vector connecting start and end of cylinder
    double proj_between = proj - start2atom;
    // This tells us if we are outside the end of the cylinder
    double excess = proj_between - d1_len;
    // Return if we are outside of the cylinder as calculated based on excess
    if( excess>low_sf.get_dmax() || -proj_between>low_sf.get_dmax() ) continue;
    // Calculate the excess swiching functions
    double edf1, eval1 = low_sf.calculate( excess, edf1 );
    double edf2, eval2 = low_sf.calculate( -proj_between, edf2 );
    // Calculate the projection on the perpendicular distance from the center of the tube
    double cm = dstart.modulo2() - proj*proj;

    // Now calculate the density in the cylinder
    if( cm>0 && cm<cylinder_sw.get_dmax2() ) {
      double dfuncr, val = cylinder_sw.calculateSqr( cm, dfuncr );
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
        dd1 = matmul(-dstart, d1_a1) - 0.5*d1;
        dd2 = matmul(-dstart, -d1_a1) - 0.5*d1;
        dd3 = d1;

        // Calculate derivatives of cross product
        Vector der( -0.5*binlength*matmul( d1_a1,dstart ) );
        dc1 = dfuncr*( 0.5*dstart + der - proj*dd1 );
        dc2 = dfuncr*( 0.5*dstart - der - proj*dd2 );
        dc3 = dfuncr*( -dstart - proj*dd3 );

        // Calculate derivatives of excess
        de1 = eval2*edf1*excess*(dd1 + 0.5*d1 ) + eval1*edf2*proj_between*(dd1 - 0.5*d1);
        de2 = eval2*edf1*excess*(dd2 - 0.5*d1 ) + eval1*edf2*proj_between*(dd2 + 0.5*d1);
        de3 = ( eval2*edf1*excess + eval1*edf2*proj_between )*dd3;
      }
      for(unsigned bin=0; bin<maxbins; ++bin) {
        bead.set( bin*binw_mat, (bin+1)*binw_mat, sigma );
        if( proj<(bin*binw_mat-bead.getCutoff()) || proj>binw_mat*(bin+1)+bead.getCutoff() ) continue;
        double der, contr=bead.calculateWithCutoff( proj, der ) / cell_volume; der /= cell_volume;
        tvals.addValue( bin, contr*val*eval1*eval2 );

        if( !doNotCalculateDerivatives() ) {
          g1derivf=contr*eval1*eval2*dc1 + val*eval1*eval2*der*dd1 + contr*val*de1;
          tvals.addDerivative( bin, 3*myvals.getTaskIndex()+0, g1derivf[0] );
          tvals.addDerivative( bin, 3*myvals.getTaskIndex()+1, g1derivf[1] );
          tvals.addDerivative( bin, 3*myvals.getTaskIndex()+2, g1derivf[2] );
          g2derivf=contr*eval1*eval2*dc2 + val*eval1*eval2*der*dd2 + contr*val*de2;
          tvals.addDerivative( bin, 3*myvals.getSecondTaskIndex()+0, g2derivf[0] );
          tvals.addDerivative( bin, 3*myvals.getSecondTaskIndex()+1, g2derivf[1] );
          tvals.addDerivative( bin, 3*myvals.getSecondTaskIndex()+2, g2derivf[2] );
          lderivf=contr*eval1*eval2*dc3 + val*eval1*eval2*der*dd3 + contr*val*de3;
          unsigned tindex = myvals.getIndices()[ i + myvals.getSplitIndex() ];
          tvals.addDerivative( bin, 3*tindex+0, lderivf[0] );
          tvals.addDerivative( bin, 3*tindex+1, lderivf[1] );
          tvals.addDerivative( bin, 3*tindex+2, lderivf[2] );
          // Virial
          vir = - Tensor( d20, g1derivf ) - Tensor( d21, g2derivf );
          unsigned nbase = 3*getNumberOfAtoms();
          tvals.addDerivative( bin, nbase+0, vir(0,0) );
          tvals.addDerivative( bin, nbase+1, vir(0,1) );
          tvals.addDerivative( bin, nbase+2, vir(0,2) );
          tvals.addDerivative( bin, nbase+3, vir(1,0) );
          tvals.addDerivative( bin, nbase+4, vir(1,1) );
          tvals.addDerivative( bin, nbase+5, vir(1,2) );
          tvals.addDerivative( bin, nbase+6, vir(2,0) );
          tvals.addDerivative( bin, nbase+7, vir(2,1) );
          tvals.addDerivative( bin, nbase+8, vir(2,2) );
        }
      }
    }
  }
  // Find maximum density
  double max = tvals.get(0); unsigned vout = 0;
  for(unsigned i=1; i<maxbins; ++i) {
    if( tvals.get(i)>max ) { max=tvals.get(i); vout=i; }
  }
  // Transform the density
  double df, tsw = threshold_switch.calculate( max, df );
  if( fabs(sw*tsw)<epsilon ) return 0;

  if( !doNotCalculateDerivatives() ) {
    Vector ader; Tensor vir; Vector ddd = tsw*dfuncl*distance;
    ader[0] = tvals.getDerivative( vout, 3*myvals.getTaskIndex()+0 );
    ader[1] = tvals.getDerivative( vout, 3*myvals.getTaskIndex()+1 );
    ader[2] = tvals.getDerivative( vout, 3*myvals.getTaskIndex()+2 );
    addAtomDerivatives( 0, sw*df*max*ader - ddd, myvals );
    ader[0] = tvals.getDerivative( vout, 3*myvals.getSecondTaskIndex()+0 );
    ader[1] = tvals.getDerivative( vout, 3*myvals.getSecondTaskIndex()+1 );
    ader[2] = tvals.getDerivative( vout, 3*myvals.getSecondTaskIndex()+2 );
    addAtomDerivatives( 1, sw*df*max*ader + ddd, myvals );
    for(unsigned i=0; i<natoms; ++i) {
      unsigned tindex = myvals.getIndices()[ i + myvals.getSplitIndex() ];
      ader[0] = tvals.getDerivative( vout, 3*tindex+0 );
      ader[1] = tvals.getDerivative( vout, 3*tindex+1 );
      ader[2] = tvals.getDerivative( vout, 3*tindex+2 );
      addThirdAtomDerivatives( i, sw*df*max*ader, myvals );
    }
    unsigned nbase = 3*getNumberOfAtoms(); Tensor vird(ddd,distance);
    vir(0,0) = sw*df*max*tvals.getDerivative( vout, nbase+0 ) - vird(0,0);
    vir(0,1) = sw*df*max*tvals.getDerivative( vout, nbase+1 ) - vird(0,1);
    vir(0,2) = sw*df*max*tvals.getDerivative( vout, nbase+2 ) - vird(0,2);
    vir(1,0) = sw*df*max*tvals.getDerivative( vout, nbase+3 ) - vird(1,0);
    vir(1,1) = sw*df*max*tvals.getDerivative( vout, nbase+4 ) - vird(1,1);
    vir(1,2) = sw*df*max*tvals.getDerivative( vout, nbase+5 ) - vird(1,2);
    vir(2,0) = sw*df*max*tvals.getDerivative( vout, nbase+6 ) - vird(2,0);
    vir(2,1) = sw*df*max*tvals.getDerivative( vout, nbase+7 ) - vird(2,1);
    vir(2,2) = sw*df*max*tvals.getDerivative( vout, nbase+8 ) - vird(2,2);
    addBoxDerivatives( vir, myvals );
  }
  return sw*tsw;
}

}
}
