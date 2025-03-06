/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Units.h"
#include "tools/Pbc.h"
#include "ActionVolume.h"
#include "tools/HistogramBead.h"
#include "VolumeShortcut.h"

//+PLUMEDOC VOLUMES CAVITY
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a box defined by the positions of four atoms.

Each of the base quantities calculated by a multicolvar can can be assigned to a particular point in three
dimensional space. For example, if we have the coordination numbers for all the atoms in the
system each coordination number can be assumed to lie on the position of the central atom.
Because each base quantity can be assigned to a particular point in space we can calculate functions of the
distribution of base quantities in a particular part of the box by using:

\f[
\overline{s}_{\tau} = \frac{ \sum_i f(s_i) w(u_i,v_i,w_i) }{ \sum_i w(u_i,v_i,w_i) }
\f]

where the sum is over the collective variables, \f$s_i\f$, each of which can be thought to be at \f$ (u_i,v_i,z_i)\f$.
The function \f$(s_i)\f$ can be any of the usual LESS_THAN, MORE_THAN, WITHIN etc that are used in all other multicolvars.
Notice that here (at variance with what is done in \ref AROUND) we have transformed from the usual \f$(x_i,y_i,z_i)\f$
position to a position in \f$ (u_i,v_i,z_i)\f$.  This is done using a rotation matrix as follows:

\f[
\left(
\begin{matrix}
 u_i \\
 v_i \\
 w_i
\end{matrix}
\right) = \mathbf{R}
\left(
\begin{matrix}
 x_i - x_o \\
 y_i - y_o \\
 z_i - z_o
\end{matrix}
\right)
\f]

where \f$\mathbf{R}\f$ is a rotation matrix that is calculated by constructing a set of three orthonormal vectors from the
reference positions specified by the user. The first of these unit vectors points from the first reference atom to the second.
The second is then the normal to the plane containing atoms 1,2 and 3 and the the third is the unit vector orthogonal to
these first two vectors.  \f$(x_o,y_o,z_o)\f$, meanwhile, specifies the position of the first reference atom.

In the previous function \f$ w(u_i,v_i,w_i) \f$ measures whether or not the system is in the subregion of interest. It
is equal to:

\f[
w(u_i,v_i,w_i) = \int_{0}^{u'} \int_{0}^{v'} \int_{0}^{w'} \textrm{d}u\textrm{d}v\textrm{d}w
   K\left( \frac{u - u_i}{\sigma} \right)K\left( \frac{v - v_i}{\sigma} \right)K\left( \frac{w - w_i}{\sigma} \right)
\f]

where \f$K\f$ is one of the kernel functions described on \ref histogrambead and \f$\sigma\f$ is a bandwidth parameter.
The vector connecting atom 1 to atom 4 is used to define the extent of the box in each of the \f$u\f$, \f$v\f$ and \f$w\f$
directions.  Essentially the vector connecting atom 1 to atom 4 is projected onto the three unit vectors
described above and the resulting projections determine the \f$u'\f$, \f$v'\f$ and \f$w'\f$ parameters in the above expression.

\par Examples

The following commands tell plumed to calculate the number of atoms in an ion channel in a protein.
The extent of the channel is calculated from the positions of atoms 1, 4, 5 and 11. The final value will be labeled cav.

\plumedfile
d1: DENSITY SPECIES=20-500
CAVITY DATA=d1 ATOMS=1,4,5,11 SIGMA=0.1 LABEL=cav
\endplumedfile

The following command tells plumed to calculate the coordination numbers (with other water molecules) for the water
molecules in the protein channel described above.  The average coordination number and the number of coordination
numbers more than 4 is then calculated.  The values of these two quantities are given the labels cav.mean and cav.morethan

\plumedfile
d1: COORDINATIONNUMBER SPECIES=20-500 R_0=0.1
CAVITY DATA=d1 ATOMS=1,4,5,11 SIGMA=0.1 MEAN MORE_THAN={RATIONAL R_0=4} LABEL=cav
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR CAVITY_CALC
/*
Calculate a vector from the input positions with elements equal to one when the positions are in a particular part of the cell and elements equal to zero otherwise

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {

class VolumeCavity {
public:
  double jacob_det;
  double len_bi, len_cross, len_perp, sigma;
  Vector bi, cross, perp;
  std::string kerneltype;
  std::vector<Vector> dlbi, dlcross, dlperp;
  std::vector<Tensor> dbi, dcross, dperp;
  static void registerKeywords( Keywords& keys );
  VolumeCavity() : jacob_det(0), len_bi(0), len_cross(0), len_perp(0), sigma(0), dlbi(4), dlcross(4), dlperp(4), dbi(3), dcross(3), dperp(3) {}
  void setupRegions( ActionVolume<VolumeCavity>* action, const Pbc& pbc, const std::vector<Vector>& positions );
  void parseInput( ActionVolume<VolumeCavity>* action );
  static void parseAtoms( ActionVolume<VolumeCavity>* action, std::vector<AtomNumber>& atoms );
  VolumeCavity& operator=( const VolumeCavity& m ) {
    jacob_det=m.jacob_det;
    len_bi=m.len_bi;
    len_cross=m.len_cross;
    len_perp=m.len_perp;
    sigma=m.sigma;
    dlbi.resize(4);
    dlcross.resize(4);
    dlperp.resize(4);
    dbi.resize(3);
    dcross.resize(3);
    dperp.resize(3);
    kerneltype=m.kerneltype;
    return *this;
  }
  static void calculateNumberInside( const VolumeInput& input, const VolumeCavity& actioninput, VolumeOutput& output );
};

typedef ActionVolume<VolumeCavity> VolCav;
PLUMED_REGISTER_ACTION(VolCav,"CAVITY_CALC")
char glob_cavity[] = "CAVITY";
typedef VolumeShortcut<glob_cavity> VolumeCavityShortcut;
PLUMED_REGISTER_ACTION(VolumeCavityShortcut,"CAVITY")

void VolumeCavity::registerKeywords( Keywords& keys ) {
  keys.setDisplayName("CAVITY");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.add("atoms","BOX","the positions of four atoms that define spatial extent of the cavity");
}

void VolumeCavity::parseInput( ActionVolume<VolumeCavity>* action ) {
  action->parse("SIGMA",sigma);
  action->parse("KERNEL",kerneltype);
  action->log.printf("  using %s kernels with a bandwidth of %d \n", kerneltype.c_str(), sigma );
}

void VolumeCavity::parseAtoms( ActionVolume<VolumeCavity>* action, std::vector<AtomNumber>& atoms ) {
  action->parseAtomList("BOX",atoms);
  if( atoms.size()!=4 ) {
    action->error("number of atoms in box should be equal to four");
  }

  action->log.printf("  boundaries for region are calculated based on positions of atoms : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    action->log.printf("%d ",atoms[i].serial() );
  }
  action->log.printf("\n");
}

void VolumeCavity::setupRegions( ActionVolume<VolumeCavity>* action, const Pbc& pbc, const std::vector<Vector>& positions ) {
  // Make some space for things
  Vector d1, d2, d3;

  // Set the position of the origin
  Vector origin=positions[0];

  // Get two vectors
  d1 = pbc.distance(origin,positions[1]);
  double d1l=d1.modulo();
  d2 = pbc.distance(origin,positions[2]);

  // Find the vector connecting the origin to the top corner of
  // the subregion
  d3 = pbc.distance(origin,positions[3]);

  // Create a set of unit vectors
  bi = d1 / d1l;
  len_bi=dotProduct( d3, bi );
  cross = crossProduct( d1, d2 );
  double crossmod=cross.modulo();
  cross = cross / crossmod;
  len_cross=dotProduct( d3, cross );
  perp = crossProduct( cross, bi );
  len_perp=dotProduct( d3, perp );

  // Calculate derivatives of box shape with respect to atoms
  double d1l3=d1l*d1l*d1l;
  dbi[0](0,0) = ( -(d1[1]*d1[1]+d1[2]*d1[2])/d1l3 );   // dx/dx
  dbi[0](0,1) = (  d1[0]*d1[1]/d1l3 );                 // dx/dy
  dbi[0](0,2) = (  d1[0]*d1[2]/d1l3 );                 // dx/dz
  dbi[0](1,0) = (  d1[1]*d1[0]/d1l3 );                 // dy/dx
  dbi[0](1,1) = ( -(d1[0]*d1[0]+d1[2]*d1[2])/d1l3 );   // dy/dy
  dbi[0](1,2) = (  d1[1]*d1[2]/d1l3 );
  dbi[0](2,0) = (  d1[2]*d1[0]/d1l3 );
  dbi[0](2,1) = (  d1[2]*d1[1]/d1l3 );
  dbi[0](2,2) = ( -(d1[1]*d1[1]+d1[0]*d1[0])/d1l3 );

  dbi[1](0,0) = ( (d1[1]*d1[1]+d1[2]*d1[2])/d1l3 );
  dbi[1](0,1) = ( -d1[0]*d1[1]/d1l3 );
  dbi[1](0,2) = ( -d1[0]*d1[2]/d1l3 );
  dbi[1](1,0) = ( -d1[1]*d1[0]/d1l3 );
  dbi[1](1,1) = ( (d1[0]*d1[0]+d1[2]*d1[2])/d1l3 );
  dbi[1](1,2) = ( -d1[1]*d1[2]/d1l3 );
  dbi[1](2,0) = ( -d1[2]*d1[0]/d1l3 );
  dbi[1](2,1) = ( -d1[2]*d1[1]/d1l3 );
  dbi[1](2,2) = ( (d1[1]*d1[1]+d1[0]*d1[0])/d1l3 );
  dbi[2].zero();

  Tensor tcderiv;
  double cmod3=crossmod*crossmod*crossmod;
  Vector ucross=crossmod*cross;
  tcderiv.setCol( 0, crossProduct( d1, Vector(-1.0,0.0,0.0) ) + crossProduct( Vector(-1.0,0.0,0.0), d2 ) );
  tcderiv.setCol( 1, crossProduct( d1, Vector(0.0,-1.0,0.0) ) + crossProduct( Vector(0.0,-1.0,0.0), d2 ) );
  tcderiv.setCol( 2, crossProduct( d1, Vector(0.0,0.0,-1.0) ) + crossProduct( Vector(0.0,0.0,-1.0), d2 ) );
  dcross[0](0,0)=( tcderiv(0,0)/crossmod - ucross[0]*(ucross[0]*tcderiv(0,0) + ucross[1]*tcderiv(1,0) + ucross[2]*tcderiv(2,0))/cmod3 );    // dx/dx
  dcross[0](0,1)=( tcderiv(0,1)/crossmod - ucross[0]*(ucross[0]*tcderiv(0,1) + ucross[1]*tcderiv(1,1) + ucross[2]*tcderiv(2,1))/cmod3 );    // dx/dy
  dcross[0](0,2)=( tcderiv(0,2)/crossmod - ucross[0]*(ucross[0]*tcderiv(0,2) + ucross[1]*tcderiv(1,2) + ucross[2]*tcderiv(2,2))/cmod3 );    // dx/dz
  dcross[0](1,0)=( tcderiv(1,0)/crossmod - ucross[1]*(ucross[0]*tcderiv(0,0) + ucross[1]*tcderiv(1,0) + ucross[2]*tcderiv(2,0))/cmod3 );    // dy/dx
  dcross[0](1,1)=( tcderiv(1,1)/crossmod - ucross[1]*(ucross[0]*tcderiv(0,1) + ucross[1]*tcderiv(1,1) + ucross[2]*tcderiv(2,1))/cmod3 );    // dy/dy
  dcross[0](1,2)=( tcderiv(1,2)/crossmod - ucross[1]*(ucross[0]*tcderiv(0,2) + ucross[1]*tcderiv(1,2) + ucross[2]*tcderiv(2,2))/cmod3 );    // dy/dz
  dcross[0](2,0)=( tcderiv(2,0)/crossmod - ucross[2]*(ucross[0]*tcderiv(0,0) + ucross[1]*tcderiv(1,0) + ucross[2]*tcderiv(2,0))/cmod3 );    // dz/dx
  dcross[0](2,1)=( tcderiv(2,1)/crossmod - ucross[2]*(ucross[0]*tcderiv(0,1) + ucross[1]*tcderiv(1,1) + ucross[2]*tcderiv(2,1))/cmod3 );    // dz/dy
  dcross[0](2,2)=( tcderiv(2,2)/crossmod - ucross[2]*(ucross[0]*tcderiv(0,2) + ucross[1]*tcderiv(1,2) + ucross[2]*tcderiv(2,2))/cmod3 );    // dz/dz

  tcderiv.setCol( 0, crossProduct( Vector(1.0,0.0,0.0), d2 ) );
  tcderiv.setCol( 1, crossProduct( Vector(0.0,1.0,0.0), d2 ) );
  tcderiv.setCol( 2, crossProduct( Vector(0.0,0.0,1.0), d2 ) );
  dcross[1](0,0)=( tcderiv(0,0)/crossmod - ucross[0]*(ucross[0]*tcderiv(0,0) + ucross[1]*tcderiv(1,0) + ucross[2]*tcderiv(2,0))/cmod3 );    // dx/dx
  dcross[1](0,1)=( tcderiv(0,1)/crossmod - ucross[0]*(ucross[0]*tcderiv(0,1) + ucross[1]*tcderiv(1,1) + ucross[2]*tcderiv(2,1))/cmod3 );    // dx/dy
  dcross[1](0,2)=( tcderiv(0,2)/crossmod - ucross[0]*(ucross[0]*tcderiv(0,2) + ucross[1]*tcderiv(1,2) + ucross[2]*tcderiv(2,2))/cmod3 );    // dx/dz
  dcross[1](1,0)=( tcderiv(1,0)/crossmod - ucross[1]*(ucross[0]*tcderiv(0,0) + ucross[1]*tcderiv(1,0) + ucross[2]*tcderiv(2,0))/cmod3 );    // dy/dx
  dcross[1](1,1)=( tcderiv(1,1)/crossmod - ucross[1]*(ucross[0]*tcderiv(0,1) + ucross[1]*tcderiv(1,1) + ucross[2]*tcderiv(2,1))/cmod3 );    // dy/dy
  dcross[1](1,2)=( tcderiv(1,2)/crossmod - ucross[1]*(ucross[0]*tcderiv(0,2) + ucross[1]*tcderiv(1,2) + ucross[2]*tcderiv(2,2))/cmod3 );    // dy/dz
  dcross[1](2,0)=( tcderiv(2,0)/crossmod - ucross[2]*(ucross[0]*tcderiv(0,0) + ucross[1]*tcderiv(1,0) + ucross[2]*tcderiv(2,0))/cmod3 );    // dz/dx
  dcross[1](2,1)=( tcderiv(2,1)/crossmod - ucross[2]*(ucross[0]*tcderiv(0,1) + ucross[1]*tcderiv(1,1) + ucross[2]*tcderiv(2,1))/cmod3 );    // dz/dy
  dcross[1](2,2)=( tcderiv(2,2)/crossmod - ucross[2]*(ucross[0]*tcderiv(0,2) + ucross[1]*tcderiv(1,2) + ucross[2]*tcderiv(2,2))/cmod3 );    // dz/dz

  tcderiv.setCol( 0, crossProduct( d1, Vector(1.0,0.0,0.0) ) );
  tcderiv.setCol( 1, crossProduct( d1, Vector(0.0,1.0,0.0) ) );
  tcderiv.setCol( 2, crossProduct( d1, Vector(0.0,0.0,1.0) ) );
  dcross[2](0,0)=( tcderiv(0,0)/crossmod - ucross[0]*(ucross[0]*tcderiv(0,0) + ucross[1]*tcderiv(1,0) + ucross[2]*tcderiv(2,0))/cmod3 );    // dx/dx
  dcross[2](0,1)=( tcderiv(0,1)/crossmod - ucross[0]*(ucross[0]*tcderiv(0,1) + ucross[1]*tcderiv(1,1) + ucross[2]*tcderiv(2,1))/cmod3 );    // dx/dy
  dcross[2](0,2)=( tcderiv(0,2)/crossmod - ucross[0]*(ucross[0]*tcderiv(0,2) + ucross[1]*tcderiv(1,2) + ucross[2]*tcderiv(2,2))/cmod3 );    // dx/dz
  dcross[2](1,0)=( tcderiv(1,0)/crossmod - ucross[1]*(ucross[0]*tcderiv(0,0) + ucross[1]*tcderiv(1,0) + ucross[2]*tcderiv(2,0))/cmod3 );    // dy/dx
  dcross[2](1,1)=( tcderiv(1,1)/crossmod - ucross[1]*(ucross[0]*tcderiv(0,1) + ucross[1]*tcderiv(1,1) + ucross[2]*tcderiv(2,1))/cmod3 );    // dy/dy
  dcross[2](1,2)=( tcderiv(1,2)/crossmod - ucross[1]*(ucross[0]*tcderiv(0,2) + ucross[1]*tcderiv(1,2) + ucross[2]*tcderiv(2,2))/cmod3 );    // dy/dz
  dcross[2](2,0)=( tcderiv(2,0)/crossmod - ucross[2]*(ucross[0]*tcderiv(0,0) + ucross[1]*tcderiv(1,0) + ucross[2]*tcderiv(2,0))/cmod3 );    // dz/dx
  dcross[2](2,1)=( tcderiv(2,1)/crossmod - ucross[2]*(ucross[0]*tcderiv(0,1) + ucross[1]*tcderiv(1,1) + ucross[2]*tcderiv(2,1))/cmod3 );    // dz/dy
  dcross[2](2,2)=( tcderiv(2,2)/crossmod - ucross[2]*(ucross[0]*tcderiv(0,2) + ucross[1]*tcderiv(1,2) + ucross[2]*tcderiv(2,2))/cmod3 );    // dz/dz

  dperp[0].setCol( 0, ( crossProduct( dcross[0].getCol(0), bi ) + crossProduct( cross, dbi[0].getCol(0) ) ) );
  dperp[0].setCol( 1, ( crossProduct( dcross[0].getCol(1), bi ) + crossProduct( cross, dbi[0].getCol(1) ) ) );
  dperp[0].setCol( 2, ( crossProduct( dcross[0].getCol(2), bi ) + crossProduct( cross, dbi[0].getCol(2) ) ) );

  dperp[1].setCol( 0, ( crossProduct( dcross[1].getCol(0), bi ) + crossProduct( cross, dbi[1].getCol(0) ) ) );
  dperp[1].setCol( 1, ( crossProduct( dcross[1].getCol(1), bi ) + crossProduct( cross, dbi[1].getCol(1) ) ) );
  dperp[1].setCol( 2, ( crossProduct( dcross[1].getCol(2), bi ) + crossProduct( cross, dbi[1].getCol(2) ) ) );

  dperp[2].setCol( 0, ( crossProduct( dcross[2].getCol(0), bi ) ) );
  dperp[2].setCol( 1, ( crossProduct( dcross[2].getCol(1), bi ) ) );
  dperp[2].setCol( 2, ( crossProduct( dcross[2].getCol(2), bi ) ) );

  // Ensure that all lengths are positive
  if( len_bi<0 ) {
    bi=-bi;
    len_bi=-len_bi;
    for(unsigned i=0; i<3; ++i) {
      dbi[i]*=-1.0;
    }
  }
  if( len_cross<0 ) {
    cross=-cross;
    len_cross=-len_cross;
    for(unsigned i=0; i<3; ++i) {
      dcross[i]*=-1.0;
    }
  }
  if( len_perp<0 ) {
    perp=-perp;
    len_perp=-len_perp;
    for(unsigned i=0; i<3; ++i) {
      dperp[i]*=-1.0;
    }
  }
  if( len_bi<=0 || len_cross<=0 || len_bi<=0 ) {
    plumed_merror("Invalid box coordinates");
  }

  // Now derivatives of lengths
  Tensor dd3( Tensor::identity() );
  dlbi[0] = matmul(d3,dbi[0]) - matmul(bi,dd3);
  dlbi[1] = matmul(d3,dbi[1]);
  dlbi[2] = matmul(d3,dbi[2]);
  dlbi[3] = matmul(bi,dd3);

  dlcross[0] = matmul(d3,dcross[0]) - matmul(cross,dd3);
  dlcross[1] = matmul(d3,dcross[1]);
  dlcross[2] = matmul(d3,dcross[2]);
  dlcross[3] = matmul(cross,dd3);

  dlperp[0] = matmul(d3,dperp[0]) - matmul(perp,dd3);
  dlperp[1] = matmul(d3,dperp[1]);
  dlperp[2] = matmul(d3,dperp[2]);
  dlperp[3] = matmul(perp,dd3);

  // Need to calculate the jacobian
  Tensor jacob;
  jacob(0,0)=bi[0];
  jacob(1,0)=bi[1];
  jacob(2,0)=bi[2];
  jacob(0,1)=cross[0];
  jacob(1,1)=cross[1];
  jacob(2,1)=cross[2];
  jacob(0,2)=perp[0];
  jacob(1,2)=perp[1];
  jacob(2,2)=perp[2];
  jacob_det = fabs( jacob.determinant() );
}

void VolumeCavity::calculateNumberInside( const VolumeInput& input, const VolumeCavity& actioninput, VolumeOutput& output ) {
  // Setup the histogram bead
  HistogramBead bead;
  bead.isNotPeriodic();
  bead.setKernelType( actioninput.kerneltype );

  // Calculate distance of atom from origin of new coordinate frame
  Vector datom=input.pbc.distance( Vector(input.refpos[0][0],input.refpos[0][1],input.refpos[0][2]), Vector(input.cpos[0],input.cpos[1],input.cpos[2]) );
  double ucontr, uder, vcontr, vder, wcontr, wder;

  // Calculate contribution from integral along bi
  bead.set( 0, actioninput.len_bi, actioninput.sigma );
  double upos=dotProduct( datom, actioninput.bi );
  ucontr=bead.calculate( upos, uder );
  double udlen=bead.uboundDerivative( upos );
  double uder2 = bead.lboundDerivative( upos ) - udlen;

  // Calculate contribution from integral along cross
  bead.set( 0, actioninput.len_cross, actioninput.sigma );
  double vpos=dotProduct( datom, actioninput.cross );
  vcontr=bead.calculate( vpos, vder );
  double vdlen=bead.uboundDerivative( vpos );
  double vder2 = bead.lboundDerivative( vpos ) - vdlen;

  // Calculate contribution from integral along perp
  bead.set( 0, actioninput.len_perp, actioninput.sigma );
  double wpos=dotProduct( datom, actioninput.perp );
  wcontr=bead.calculate( wpos, wder );
  double wdlen=bead.uboundDerivative( wpos );
  double wder2 = bead.lboundDerivative( wpos ) - wdlen;

  Vector dfd;
  dfd[0]=uder*vcontr*wcontr;
  dfd[1]=ucontr*vder*wcontr;
  dfd[2]=ucontr*vcontr*wder;
  output.derivatives[0] = (dfd[0]*actioninput.bi[0]+dfd[1]*actioninput.cross[0]+dfd[2]*actioninput.perp[0]);
  output.derivatives[1] = (dfd[0]*actioninput.bi[1]+dfd[1]*actioninput.cross[1]+dfd[2]*actioninput.perp[1]);
  output.derivatives[2] = (dfd[0]*actioninput.bi[2]+dfd[1]*actioninput.cross[2]+dfd[2]*actioninput.perp[2]);
  output.values[0] = ucontr*vcontr*wcontr*actioninput.jacob_det;

  // Add reference atom derivatives
  dfd[0]=uder2*vcontr*wcontr;
  dfd[1]=ucontr*vder2*wcontr;
  dfd[2]=ucontr*vcontr*wder2;
  Vector dfld;
  dfld[0]=udlen*vcontr*wcontr;
  dfld[1]=ucontr*vdlen*wcontr;
  dfld[2]=ucontr*vcontr*wdlen;
  output.refders[0] = dfd[0]*matmul(datom,actioninput.dbi[0]) + dfd[1]*matmul(datom,actioninput.dcross[0]) + dfd[2]*matmul(datom,actioninput.dperp[0]) +
                      dfld[0]*actioninput.dlbi[0] + dfld[1]*actioninput.dlcross[0] + dfld[2]*actioninput.dlperp[0] - Vector(output.derivatives[0],output.derivatives[1],output.derivatives[2]);
  output.refders[1] = dfd[0]*matmul(datom,actioninput.dbi[1]) + dfd[1]*matmul(datom,actioninput.dcross[1]) + dfd[2]*matmul(datom,actioninput.dperp[1]) +
                      dfld[0]*actioninput.dlbi[1] + dfld[1]*actioninput.dlcross[1] + dfld[2]*actioninput.dlperp[1];
  output.refders[2] = dfd[0]*matmul(datom,actioninput.dbi[2]) + dfd[1]*matmul(datom,actioninput.dcross[2]) + dfd[2]*matmul(datom,actioninput.dperp[2]) +
                      dfld[0]*actioninput.dlbi[2] + dfld[1]*actioninput.dlcross[2] + dfld[2]*actioninput.dlperp[2];
  output.refders[3] = dfld[0]*actioninput.dlbi[3] + dfld[1]*actioninput.dlcross[3] + dfld[2]*actioninput.dlperp[3];

  Tensor vir;
  vir=-Tensor( Vector(input.cpos[0],input.cpos[1],input.cpos[2]), Vector(output.derivatives[0],output.derivatives[1],output.derivatives[2]) );
  for(unsigned i=0; i<4; ++i) {
    vir -= Tensor( Vector(input.refpos[i][0],input.refpos[i][1],input.refpos[i][2]), Vector(output.refders[i][0],output.refders[i][1],output.refders[i][2]) );
  }
  output.virial.set( 0, vir );
}

}
}
