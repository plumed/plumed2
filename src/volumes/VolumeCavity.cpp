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
#include "VolumeShortcut.h"

//+PLUMEDOC VOLUMES CAVITY
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a box defined by the positions of four atoms.

This action can be used similarly to the way [AROUND](AROUND.md) is used.  Like [AROUND](AROUND.md) this action returns a vector
with elements that tell you whether an input atom is within a part of the box that is of particular interest or not. However, for this action
the elements of this vector are calculated using:

$$
w(u_i,v_i,w_i) = \int_{0}^{u'} \int_{0}^{v'} \int_{0}^{w'} \textrm{d}u\textrm{d}v\textrm{d}w
   K\left( \frac{u - u_i}{\sigma} \right)K\left( \frac{v - v_i}{\sigma} \right)K\left( \frac{w - w_i}{\sigma} \right)
$$

with $u_i,v_i,w_i$ being calculated from the positon of the $i$th atom, $(x_i,y_i,z_i)$, by performing the following transformation.

$$
\left(
\begin{matrix}
 u_i \\
 v_i \\
 w_i
\end{matrix}
\right) = R
\left(
\begin{matrix}
 x_i - x_o \\
 y_i - y_o \\
 z_i - z_o
\end{matrix}
\right)
$$

In this expression $R$ is a rotation matrix that is calculated by constructing a set of three orthonormal vectors from the
reference positions specified by the user. The first of these unit vectors points from the first reference atom to the second.
The second is then the normal to the plane containing atoms 1,2 and 3 and the the third is the unit vector orthogonal to
these first two vectors.  $(x_o,y_o,z_o)$, meanwhile, specifies the position of the first reference atom.

In the first expression above $K$ is one of the kernel functions described in the documentation for the action [BETWEEN](BETWEEN.md)
and $\sigma$ is a bandwidth parameter.  Furthermore, The vector connecting atom 1 to atom 4 is used to define the extent of the box in
each of the $u$, $v$ and $w$ directions.  This vector connecting atom 1 to atom 4 is projected onto the three unit vectors
described above and the resulting projections determine the $u'$, $v'$ and $w'$ parameters in the above expression.

The following commands illustrate how this works in practise.  We are using PLUMED here to calculate the number of atoms from the group specified using the ATOMS keyword below are
in an ion channel in a protein.  The extent of the channel is calculated from the positions of atoms 1, 4, 5 and 11.

```plumed
cav: CAVITY ATOMS=20-500 BOX=1,4,5,11 SIGMA=0.1
s: SUM ARG=cav PERIODIC=NO
PRINT ARG=s FILE=colvar
```

By contrst the following command tells plumed to calculate the coordination numbers (with other water molecules) for the water
molecules in the protein channel described above. The average coordination number and the number of coordination
numbers more than 4 is then calculated for those molecules that are in the region of interest.

```plumed
# Calculate the atoms that are in the cavity
cav: CAVITY ATOMS=20-500 BOX=1,4,5,11 SIGMA=0.1
# Calculate the coordination numbers of all the atoms
d1: COORDINATIONNUMBER SPECIES=20-500 SWITCH={RATIONAL R_0=0.1}
# Do atoms have a coordination number that is greater than 4
fd1: MORE_THAN ARG=d1 SWITCH={RATIONAL R_0=4}
# Calculate the mean coordination number in the channel
nnn: CUSTOM ARG=cav,d1 FUNC=x*y PERIODIC=NO
numer: SUM ARG=nnn PERIODIC=NO
denom: SUM ARG=cav PERIODIC=NO
mean: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# Calculate the number of atoms that are in the channel and that have a coordination number that is greater than 4
sss: CUSTOM ARG=fd1,cav FUNC=x*y PERIODIC=NO
mt: SUM ARG=sss PERIODIC=NO
# And print these two quantities to a file
PRINT ARG=mean,mt FILE=colvar
```

As with [AROUND](AROUND.md) earlier version of PLUMED used a different syntax for doing these types of calculations, which can
still be used with this new version of the command.  However, we strongly recommend using the newer syntax.

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

class VolumeCavity : public ActionVolume {
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
  explicit VolumeCavity(const ActionOptions& ao);
  ~VolumeCavity();
  void setupRegions() override;
  void update() override;
  double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const override;
};

PLUMED_REGISTER_ACTION(VolumeCavity,"CAVITY_CALC")
char glob_cavity[] = "CAVITY";
typedef VolumeShortcut<glob_cavity> VolumeCavityShortcut;
PLUMED_REGISTER_ACTION(VolumeCavityShortcut,"CAVITY")

void VolumeCavity::registerKeywords( Keywords& keys ) {
  ActionVolume::registerKeywords( keys );
  keys.setDisplayName("CAVITY");
  keys.add("atoms","BOX","the positions of four atoms that define spatial extent of the cavity");
  keys.addFlag("PRINT_BOX",false,"write out the positions of the corners of the box to an xyz file");
  keys.add("optional","FILE","the file on which to write out the box coordinates");
  keys.add("optional","UNITS","( default=nm ) the units in which to write out the corners of the box");
}

VolumeCavity::VolumeCavity(const ActionOptions& ao):
  Action(ao),
  ActionVolume(ao),
  boxout(false),
  lenunit(1.0),
  dlbi(4),
  dlcross(4),
  dlperp(4),
  dbi(3),
  dcross(3),
  dperp(3) {
  std::vector<AtomNumber> atoms;
  parseAtomList("BOX",atoms);
  if( atoms.size()!=4 ) {
    error("number of atoms in box should be equal to four");
  }

  log.printf("  boundaries for region are calculated based on positions of atoms : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    log.printf("%d ",atoms[i].serial() );
  }
  log.printf("\n");
  requestAtoms( atoms );

  boxout=false;
  parseFlag("PRINT_BOX",boxout);
  if(boxout) {
    std::string boxfname;
    parse("FILE",boxfname);
    if(boxfname.length()==0) {
      error("no name for box file specified");
    }
    std::string unitname;
    parse("UNITS",unitname);
    if ( unitname.length()>0 ) {
      Units u;
      u.setLength(unitname);
      lenunit=getUnits().getLength()/u.getLength();
    } else {
      unitname="nm";
    }
    boxfile.link(*this);
    boxfile.open( boxfname );
    log.printf("  printing box coordinates on file named %s in %s \n",boxfname.c_str(), unitname.c_str() );
  }

  checkRead();
}

VolumeCavity::~VolumeCavity() {
}

void VolumeCavity::setupRegions() {
  // Make some space for things
  Vector d1, d2, d3;

  // Retrieve the sigma value
  sigma=getSigma();
  // Set the position of the origin
  origin=getPosition(0);

  // Get two vectors
  d1 = pbcDistance(origin,getPosition(1));
  double d1l=d1.modulo();
  d2 = pbcDistance(origin,getPosition(2));

  // Find the vector connecting the origin to the top corner of
  // the subregion
  d3 = pbcDistance(origin,getPosition(3));

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

void VolumeCavity::update() {
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

double VolumeCavity::calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& rderiv ) const {
  // Setup the histogram bead
  HistogramBead bead;
  bead.isNotPeriodic();
  bead.setKernelType( getKernelType() );

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

  Vector dfd;
  dfd[0]=uder*vcontr*wcontr;
  dfd[1]=ucontr*vder*wcontr;
  dfd[2]=ucontr*vcontr*wder;
  derivatives[0] = (dfd[0]*bi[0]+dfd[1]*cross[0]+dfd[2]*perp[0]);
  derivatives[1] = (dfd[0]*bi[1]+dfd[1]*cross[1]+dfd[2]*perp[1]);
  derivatives[2] = (dfd[0]*bi[2]+dfd[1]*cross[2]+dfd[2]*perp[2]);
  double tot = ucontr*vcontr*wcontr*jacob_det;

  // Add reference atom derivatives
  dfd[0]=uder2*vcontr*wcontr;
  dfd[1]=ucontr*vder2*wcontr;
  dfd[2]=ucontr*vcontr*wder2;
  Vector dfld;
  dfld[0]=udlen*vcontr*wcontr;
  dfld[1]=ucontr*vdlen*wcontr;
  dfld[2]=ucontr*vcontr*wdlen;
  rderiv[0] = dfd[0]*matmul(datom,dbi[0]) + dfd[1]*matmul(datom,dcross[0]) + dfd[2]*matmul(datom,dperp[0]) +
              dfld[0]*dlbi[0] + dfld[1]*dlcross[0] + dfld[2]*dlperp[0] - derivatives;
  rderiv[1] = dfd[0]*matmul(datom,dbi[1]) + dfd[1]*matmul(datom,dcross[1]) + dfd[2]*matmul(datom,dperp[1]) +
              dfld[0]*dlbi[1] + dfld[1]*dlcross[1] + dfld[2]*dlperp[1];
  rderiv[2] = dfd[0]*matmul(datom,dbi[2]) + dfd[1]*matmul(datom,dcross[2]) + dfd[2]*matmul(datom,dperp[2]) +
              dfld[0]*dlbi[2] + dfld[1]*dlcross[2] + dfld[2]*dlperp[2];
  rderiv[3] = dfld[0]*dlbi[3] + dfld[1]*dlcross[3] + dfld[2]*dlperp[3];

  vir.zero();
  vir-=Tensor( cpos,derivatives );
  for(unsigned i=0; i<4; ++i) {
    vir -= Tensor( getPosition(i), rderiv[i] );
  }

  return tot;
}

}
}
