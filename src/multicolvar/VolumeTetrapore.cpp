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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include "tools/Pbc.h"
#include "ActionVolume.h"

//+PLUMEDOC VOLUMES TETRAHEDRALPORE
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms lie that lie in a box defined by the positions of four atoms at the corners of a tetrahedron.

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
reference positions specified by the user.  Initially unit vectors are found by calculating the bisector, \f$\mathbf{b}\f$, and
cross product, \f$\mathbf{c}\f$, of the vectors connecting atoms 1 and 2.  A third unit vector, \f$\mathbf{p}\f$ is then found by taking the cross
product between the cross product calculated during the first step, \f$\mathbf{c}\f$ and the bisector, \f$\mathbf{b}\f$.  From this
second cross product \f$\mathbf{p}\f$ and the bisector \f$\mathbf{b}\f$ two new vectors are calculated using:

\f[
v_1 = \cos\left(\frac{\pi}{4}\right)\mathbf{b} + \sin\left(\frac{\pi}{4}\right)\mathbf{p} \qquad \textrm{and} \qquad
v_2 = \cos\left(\frac{\pi}{4}\right)\mathbf{b} - \sin\left(\frac{\pi}{4}\right)\mathbf{p}
\f]

In the previous function \f$ w(u_i,v_i,w_i) \f$ measures whether or not the system is in the subregion of interest. It
is equal to:

\f[
w(u_i,v_i,w_i) = \int_{0}^{u'} \int_{0}^{v'} \int_{0}^{w'} \textrm{d}u\textrm{d}v\textrm{d}w
   K\left( \frac{u - u_i}{\sigma} \right)K\left( \frac{v - v_i}{\sigma} \right)K\left( \frac{w - w_i}{\sigma} \right)
\f]

where \f$K\f$ is one of the kernel functions described on \ref histogrambead and \f$\sigma\f$ is a bandwidth parameter.
The values of \f$u'\f$ and \f$v'\f$ are found by finding the projections of the vectors connecting atoms 1 and 2 and 1
and 3 \f$v_1\f$ and \f$v_2\f$.  This gives four projections: the largest two projections are used in the remainder of
the calculations.  \f$w'\f$ is calculated by taking the projection of the vector connecting atoms 1 and 4 on the vector
\f$\mathbf{c}\f$.  Notice that the manner by which this box is constructed differs from the way this is done in \ref CAVITY.
This is in fact the only point of difference between these two actions.

\par Examples

The following commands tell plumed to calculate the number of atom inside a tetrahedral cavity.  The extent of the tetrahedral
cavity is calculated from the positions of atoms 1, 4, 5, and 11,  The final value will be labeled cav.

\plumedfile
d1: DENSITY SPECIES=20-500
TETRAHEDRALPORE DATA=d1 ATOMS=1,4,5,11 SIGMA=0.1 LABEL=cav
\endplumedfile

The following command tells plumed to calculate the coordination numbers (with other water molecules) for the water
molecules in the tetrahedral cavity described above.  The average coordination number and the number of coordination
numbers more than 4 is then calculated.  The values of these two quantities are given the labels cav.mean and cav.morethan

\plumedfile
d1: COORDINATIONNUMBER SPECIES=20-500 R_0=0.1
CAVITY DATA=d1 ATOMS=1,4,5,11 SIGMA=0.1 MEAN MORE_THAN={RATIONAL R_0=4} LABEL=cav
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class VolumeTetrapore : public ActionVolume {
private:
  bool boxout;
  OFile boxfile;
  double lenunit;
  double jacob_det;
  double len_bi, len_cross, len_perp, sigma;
  Vector origin, bi, cross, perp;
  std::vector<Vector> dlbi, dlcross, dlperp;
  std::vector<Tensor> dbi, dcross, dperp;
public:
  static void registerKeywords( Keywords& keys );
  explicit VolumeTetrapore(const ActionOptions& ao);
  ~VolumeTetrapore();
  void setupRegions() override;
  void update() override;
  double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const override;
};

PLUMED_REGISTER_ACTION(VolumeTetrapore,"TETRAHEDRALPORE")

void VolumeTetrapore::registerKeywords( Keywords& keys ) {
  ActionVolume::registerKeywords( keys );
  keys.add("atoms","ATOMS","the positions of four atoms that define spatial extent of the cavity");
  keys.addFlag("PRINT_BOX",false,"write out the positions of the corners of the box to an xyz file");
  keys.add("optional","FILE","the file on which to write out the box coordinates");
  keys.add("optional","UNITS","( default=nm ) the units in which to write out the corners of the box");
}

VolumeTetrapore::VolumeTetrapore(const ActionOptions& ao):
  Action(ao),
  ActionVolume(ao),
  boxout(false),
  lenunit(1.0),
  dlbi(4),
  dlcross(4),
  dlperp(4),
  dbi(3),
  dcross(3),
  dperp(3)
{
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if( atoms.size()!=4 ) error("number of atoms should be equal to four");

  log.printf("  boundaries for region are calculated based on positions of atoms : ");
  for(unsigned i=0; i<atoms.size(); ++i) log.printf("%d ",atoms[i].serial() );
  log.printf("\n");

  boxout=false; parseFlag("PRINT_BOX",boxout);
  if(boxout) {
    std::string boxfname; parse("FILE",boxfname);
    if(boxfname.length()==0) error("no name for box file specified");
    std::string unitname; parse("UNITS",unitname);
    if ( unitname.length()>0 ) {
      Units u; u.setLength(unitname);
      lenunit=plumed.getAtoms().getUnits().getLength()/u.getLength();
    } else {
      unitname="nm";
    }
    boxfile.link(*this);
    boxfile.open( boxfname.c_str() );
    log.printf("  printing box coordinates on file named %s in %s \n",boxfname.c_str(), unitname.c_str() );
  }

  checkRead();
  requestAtoms(atoms);
  // We have to readd the dependency because requestAtoms removes it
  addDependency( getPntrToMultiColvar() );
}

VolumeTetrapore::~VolumeTetrapore() {
}

void VolumeTetrapore::setupRegions() {
  // Make some space for things
  Vector d1, d2, d3;

  // Retrieve the sigma value
  sigma=getSigma();
  // Set the position of the origin
  origin=getPosition(0);

  // Get two vectors
  d1 = pbcDistance(origin,getPosition(1));
  d2 = pbcDistance(origin,getPosition(2));

  // Find the vector connecting the origin to the top corner of
  // the subregion
  d3 = pbcDistance(origin,getPosition(3));

  // Create a set of unit vectors
  Vector bisector = d1 + d2; double bmod=bisector.modulo(); bisector=bisector/bmod;

  // bi = d1 / d1l; len_bi=dotProduct( d3, bi );
  cross = crossProduct( d1, d2 ); double crossmod=cross.modulo();
  cross = cross / crossmod; len_cross=dotProduct( d3, cross );
  Vector truep = crossProduct( cross, bisector );

  // These are our true vectors 45 degrees from bisector
  bi = cos(pi/4.0)*bisector + sin(pi/4.0)*truep;
  perp = cos(pi/4.0)*bisector - sin(pi/4.0)*truep;

  // And the lengths of the various parts average distance to opposite corners of tetetrahedron
  len_bi = dotProduct( d1, bi ); double len_bi2 = dotProduct( d2, bi ); unsigned lbi=1;
  if( len_bi2>len_bi ) { len_bi=len_bi2; lbi=2; }
  len_perp = dotProduct( d1, perp ); double len_perp2 = dotProduct( d2, perp ); unsigned lpi=1;
  if( len_perp2>len_perp ) { len_perp=len_perp2; lpi=2; }
  plumed_assert( lbi!=lpi );

  Tensor tcderiv; double cmod3=crossmod*crossmod*crossmod; Vector ucross=crossmod*cross;
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

  std::vector<Tensor> dbisector(3);
  double bmod3=bmod*bmod*bmod; Vector ubisector=bmod*bisector;
  dbisector[0](0,0)= -2.0/bmod + 2*ubisector[0]*ubisector[0]/bmod3;
  dbisector[0](0,1)= 2*ubisector[0]*ubisector[1]/bmod3;
  dbisector[0](0,2)= 2*ubisector[0]*ubisector[2]/bmod3;
  dbisector[0](1,0)= 2*ubisector[1]*ubisector[0]/bmod3;
  dbisector[0](1,1)= -2.0/bmod + 2*ubisector[1]*ubisector[1]/bmod3;
  dbisector[0](1,2)= 2*ubisector[1]*ubisector[2]/bmod3;
  dbisector[0](2,0)= 2*ubisector[2]*ubisector[0]/bmod3;
  dbisector[0](2,1)= 2*ubisector[2]*ubisector[1]/bmod3;
  dbisector[0](2,2)= -2.0/bmod + 2*ubisector[2]*ubisector[2]/bmod3;

  dbisector[1](0,0)= 1.0/bmod - ubisector[0]*ubisector[0]/bmod3;
  dbisector[1](0,1)= -ubisector[0]*ubisector[1]/bmod3;
  dbisector[1](0,2)= -ubisector[0]*ubisector[2]/bmod3;
  dbisector[1](1,0)= -ubisector[1]*ubisector[0]/bmod3;
  dbisector[1](1,1)= 1.0/bmod - ubisector[1]*ubisector[1]/bmod3;
  dbisector[1](1,2)= -ubisector[1]*ubisector[2]/bmod3;
  dbisector[1](2,0)= -ubisector[2]*ubisector[0]/bmod3;
  dbisector[1](2,1)= -ubisector[2]*ubisector[1]/bmod3;
  dbisector[1](2,2)=1.0/bmod - ubisector[2]*ubisector[2]/bmod3;

  dbisector[2](0,0)=1.0/bmod - ubisector[0]*ubisector[0]/bmod3;
  dbisector[2](0,1)= -ubisector[0]*ubisector[1]/bmod3;
  dbisector[2](0,2)= -ubisector[0]*ubisector[2]/bmod3;
  dbisector[2](1,0)= -ubisector[1]*ubisector[0]/bmod3;
  dbisector[2](1,1)=1.0/bmod - ubisector[1]*ubisector[1]/bmod3;
  dbisector[2](1,2)= -ubisector[1]*ubisector[2]/bmod3;
  dbisector[2](2,0)= -ubisector[2]*ubisector[0]/bmod3;
  dbisector[2](2,1)= -ubisector[2]*ubisector[1]/bmod3;
  dbisector[2](2,2)=1.0/bmod - ubisector[2]*ubisector[2]/bmod3;

  std::vector<Tensor> dtruep(3);
  dtruep[0].setCol( 0, ( crossProduct( dcross[0].getCol(0), bisector ) + crossProduct( cross, dbisector[0].getCol(0) ) ) );
  dtruep[0].setCol( 1, ( crossProduct( dcross[0].getCol(1), bisector ) + crossProduct( cross, dbisector[0].getCol(1) ) ) );
  dtruep[0].setCol( 2, ( crossProduct( dcross[0].getCol(2), bisector ) + crossProduct( cross, dbisector[0].getCol(2) ) ) );

  dtruep[1].setCol( 0, ( crossProduct( dcross[1].getCol(0), bisector ) + crossProduct( cross, dbisector[1].getCol(0) ) ) );
  dtruep[1].setCol( 1, ( crossProduct( dcross[1].getCol(1), bisector ) + crossProduct( cross, dbisector[1].getCol(1) ) ) );
  dtruep[1].setCol( 2, ( crossProduct( dcross[1].getCol(2), bisector ) + crossProduct( cross, dbisector[1].getCol(2) ) ) );

  dtruep[2].setCol( 0, ( crossProduct( dcross[2].getCol(0), bisector ) + crossProduct( cross, dbisector[2].getCol(0) ) ) );
  dtruep[2].setCol( 1, ( crossProduct( dcross[2].getCol(1), bisector ) + crossProduct( cross, dbisector[2].getCol(1) ) ) );
  dtruep[2].setCol( 2, ( crossProduct( dcross[2].getCol(2), bisector ) + crossProduct( cross, dbisector[2].getCol(2) ) ) );

  // Now convert these to the derivatives of the true axis
  for(unsigned i=0; i<3; ++i) {
    dbi[i] = cos(pi/4.0)*dbisector[i] + sin(pi/4.0)*dtruep[i];
    dperp[i] = cos(pi/4.0)*dbisector[i] - sin(pi/4.0)*dtruep[i];
  }

  // Ensure that all lengths are positive
  if( len_bi<0 ) {
    bi=-bi; len_bi=-len_bi;
    for(unsigned i=0; i<3; ++i) dbi[i]*=-1.0;
  }
  if( len_cross<0 ) {
    cross=-cross; len_cross=-len_cross;
    for(unsigned i=0; i<3; ++i) dcross[i]*=-1.0;
  }
  if( len_perp<0 ) {
    perp=-perp; len_perp=-len_perp;
    for(unsigned i=0; i<3; ++i) dperp[i]*=-1.0;
  }
  if( len_bi<=0 || len_cross<=0 || len_bi<=0 ) plumed_merror("Invalid box coordinates");

  // Now derivatives of lengths
  Tensor dd3( Tensor::identity() ); Vector ddb2=d1; if( lbi==2 ) ddb2=d2;
  dlbi[1].zero(); dlbi[2].zero(); dlbi[3].zero();
  dlbi[0] = matmul(ddb2,dbi[0]) - matmul(bi,dd3);
  dlbi[lbi] = matmul(ddb2,dbi[lbi]) + matmul(bi,dd3);  // Derivative wrt d1

  dlcross[0] = matmul(d3,dcross[0]) - matmul(cross,dd3);
  dlcross[1] = matmul(d3,dcross[1]);
  dlcross[2] = matmul(d3,dcross[2]);
  dlcross[3] = matmul(cross,dd3);

  ddb2=d1; if( lpi==2 ) ddb2=d2;
  dlperp[1].zero(); dlperp[2].zero(); dlperp[3].zero();
  dlperp[0] = matmul(ddb2,dperp[0]) - matmul( perp, dd3 );
  dlperp[lpi] = matmul(ddb2,dperp[lpi]) + matmul(perp, dd3);

  // Need to calculate the jacobian
  Tensor jacob;
  jacob(0,0)=bi[0]; jacob(1,0)=bi[1]; jacob(2,0)=bi[2];
  jacob(0,1)=cross[0]; jacob(1,1)=cross[1]; jacob(2,1)=cross[2];
  jacob(0,2)=perp[0]; jacob(1,2)=perp[1]; jacob(2,2)=perp[2];
  jacob_det = fabs( jacob.determinant() );
}

void VolumeTetrapore::update() {
  if(boxout) {
    boxfile.printf("%d\n",8);
    const Tensor & t(getPbc().getBox());
    if(getPbc().isOrthorombic()) {
      boxfile.printf(" %f %f %f\n",lenunit*t(0,0),lenunit*t(1,1),lenunit*t(2,2));
    } else {
      boxfile.printf(" %f %f %f %f %f %f %f %f %f\n",
                     lenunit*t(0,0),lenunit*t(0,1),lenunit*t(0,2),
                     lenunit*t(1,0),lenunit*t(1,1),lenunit*t(1,2),
                     lenunit*t(2,0),lenunit*t(2,1),lenunit*t(2,2)
                    );
    }
    boxfile.printf("AR %f %f %f \n",lenunit*origin[0],lenunit*origin[1],lenunit*origin[2]);
    Vector ut, vt, wt;
    ut = origin + len_bi*bi;
    vt = origin + len_cross*cross;
    wt = origin + len_perp*perp;
    boxfile.printf("AR %f %f %f \n",lenunit*(ut[0]), lenunit*(ut[1]), lenunit*(ut[2]) );
    boxfile.printf("AR %f %f %f \n",lenunit*(vt[0]), lenunit*(vt[1]), lenunit*(vt[2]) );
    boxfile.printf("AR %f %f %f \n",lenunit*(wt[0]), lenunit*(wt[1]), lenunit*(wt[2]) );
    boxfile.printf("AR %f %f %f \n",lenunit*(vt[0]+len_bi*bi[0]),
                   lenunit*(vt[1]+len_bi*bi[1]),
                   lenunit*(vt[2]+len_bi*bi[2]) );
    boxfile.printf("AR %f %f %f \n",lenunit*(ut[0]+len_perp*perp[0]),
                   lenunit*(ut[1]+len_perp*perp[1]),
                   lenunit*(ut[2]+len_perp*perp[2]) );
    boxfile.printf("AR %f %f %f \n",lenunit*(vt[0]+len_perp*perp[0]),
                   lenunit*(vt[1]+len_perp*perp[1]),
                   lenunit*(vt[2]+len_perp*perp[2]) );
    boxfile.printf("AR %f %f %f \n",lenunit*(vt[0]+len_perp*perp[0]+len_bi*bi[0]),
                   lenunit*(vt[1]+len_perp*perp[1]+len_bi*bi[1]),
                   lenunit*(vt[2]+len_perp*perp[2]+len_bi*bi[2]) );
  }
}

double VolumeTetrapore::calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& rderiv ) const {
  // Setup the histogram bead
  HistogramBead bead; bead.isNotPeriodic(); bead.setKernelType( getKernelType() );

  // Calculate distance of atom from origin of new coordinate frame
  Vector datom=pbcDistance( origin, cpos );
  double ucontr, uder, vcontr, vder, wcontr, wder;

  // Calculate contribution from integral along bi
  bead.set( 0, len_bi, sigma );
  double upos=dotProduct( datom, bi );
  ucontr=bead.calculate( upos, uder );
  double udlen=bead.uboundDerivative( upos );
  double uder2 = bead.lboundDerivative( upos ) - udlen;

  // Calculate contribution from integral along cross
  bead.set( 0, len_cross, sigma );
  double vpos=dotProduct( datom, cross );
  vcontr=bead.calculate( vpos, vder );
  double vdlen=bead.uboundDerivative( vpos );
  double vder2 = bead.lboundDerivative( vpos ) - vdlen;

  // Calculate contribution from integral along perp
  bead.set( 0, len_perp, sigma );
  double wpos=dotProduct( datom, perp );
  wcontr=bead.calculate( wpos, wder );
  double wdlen=bead.uboundDerivative( wpos );
  double wder2 = bead.lboundDerivative( wpos ) - wdlen;

  Vector dfd; dfd[0]=uder*vcontr*wcontr; dfd[1]=ucontr*vder*wcontr; dfd[2]=ucontr*vcontr*wder;
  derivatives[0] = (dfd[0]*bi[0]+dfd[1]*cross[0]+dfd[2]*perp[0]);
  derivatives[1] = (dfd[0]*bi[1]+dfd[1]*cross[1]+dfd[2]*perp[1]);
  derivatives[2] = (dfd[0]*bi[2]+dfd[1]*cross[2]+dfd[2]*perp[2]);
  double tot = ucontr*vcontr*wcontr*jacob_det;

  // Add reference atom derivatives
  dfd[0]=uder2*vcontr*wcontr; dfd[1]=ucontr*vder2*wcontr; dfd[2]=ucontr*vcontr*wder2;
  Vector dfld; dfld[0]=udlen*vcontr*wcontr; dfld[1]=ucontr*vdlen*wcontr; dfld[2]=ucontr*vcontr*wdlen;
  rderiv[0] = dfd[0]*matmul(datom,dbi[0]) + dfd[1]*matmul(datom,dcross[0]) + dfd[2]*matmul(datom,dperp[0]) +
              dfld[0]*dlbi[0] + dfld[1]*dlcross[0] + dfld[2]*dlperp[0] - derivatives;
  rderiv[1] = dfd[0]*matmul(datom,dbi[1]) + dfd[1]*matmul(datom,dcross[1]) + dfd[2]*matmul(datom,dperp[1]) +
              dfld[0]*dlbi[1] + dfld[1]*dlcross[1] + dfld[2]*dlperp[1];
  rderiv[2] = dfd[0]*matmul(datom,dbi[2]) + dfd[1]*matmul(datom,dcross[2]) + dfd[2]*matmul(datom,dperp[2]) +
              dfld[0]*dlbi[2] + dfld[1]*dlcross[2] + dfld[2]*dlperp[2];
  rderiv[3] = dfld[0]*dlbi[3] + dfld[1]*dlcross[3] + dfld[2]*dlperp[3];

  vir.zero(); vir-=Tensor( cpos,derivatives );
  for(unsigned i=0; i<4; ++i) {
    vir -= Tensor( getPosition(i), rderiv[i] );
  }

  return tot;
}

}
}
