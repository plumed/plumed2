/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) crystdistrib 2023-2023 The code team
   (see the PEOPLE-crystdistrib file at the root of this folder for a list of names)

   This file is part of crystdistrib code module.

   The crystdistrib code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The crystdistrib code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the crystdistrib code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/Colvar.h"
#include "colvar/ColvarShortcut.h"
#include "colvar/MultiColvarTemplate.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"

namespace PLMD {
namespace crystdistrib {

//+PLUMEDOC COLVAR QUATERNION
/*
Calculate quaternions for molecules.

The reference frame for the molecule is defined using the positions of three user selected atom.  From the positions of these atoms, 
\f$\mathbf{x}_1\f$, \f$\mathbf{x}_2\f$ and \f$\mathbf{x}_3\f$, we define the vectors of the reference frame as:

\f[
\begin{aligned}
\mathbf{x} & = \mathbf{x}_2 - \mathbf{x}_1 \\
\mathbf{y} & = (\mathbf{x}_2 - \mathbf{x}_1) \times (\mathbf{x}_3 - \mathbf{x}_1) \\
\mathbf{z} & \mathbf{x} \times \mathbf{y} 
\f] 

\par Examples

This calculates the quaternions for a molecule with 10 atoms

\plumedfile
q1: QUATERNION ATOMS1=1,2,3
PRINT ARG=q1.w,q1.i,q1.j,q1.k FILE=colvar 
\endplumedfile

This calculate the quaternions for two molecules with 10 atoms

\plumedfile
q1: QUATERNION ATOMS1=1,2,3 ATOMS=4,5,6
PRINT ARG=q1.w,q1.i,q1.j,q1.k FILE=colvar 
\endplumedfile

*/
//+ENDPLUMEDOC

//simple hamilton/quaternion product, which is just expanding the two quats as expected, then applying the rules for i j k
//gareth - let me know if there's a way you'd prefer this be done
//passing 12 references might be a bit silly
//void QuatProduct(double &a1, double &b1, double &c1, double &d1, double &a2, double &b2, double &c2, double &d2, double &w, double &i, double &j, double &k) {
//  
//  w = a1*a2 - b1*b2 - c1*c2 - d1*d2;
//  i = a1*b2 + b1*a2 + c1*d2 - d1*c2;
//  j = a1*c2 - b1*d2 + c1*a2 + d1*b2;
//  k = a1*d2 + b1*c2 - c1*b2 + d1*a2;
//}

//
//
//

class Quaternion : public Colvar {
private:
  bool pbc;
  std::vector<double> value, masses, charges;
  std::vector<std::vector<Vector> > derivs;
  std::vector<Tensor> virial;
public:
  static void registerKeywords( Keywords& keys );
  explicit Quaternion(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                           const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                           std::vector<Tensor>& virial, const ActionAtomistic* aa );
};

typedef colvar::ColvarShortcut<Quaternion> QuaternionShortcut;
PLUMED_REGISTER_ACTION(QuaternionShortcut,"QUATERNION")
PLUMED_REGISTER_ACTION(Quaternion,"QUATERNION_SCALAR")
typedef colvar::MultiColvarTemplate<Quaternion> QuaternionMulti;
PLUMED_REGISTER_ACTION(QuaternionMulti,"QUATERNION_VECTOR")

void Quaternion::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the three atom that we are using to calculate the quaternion");
  //for working calculations, ill be using qw, qx, qy, qz, so I don't confuse myself too much 
  keys.addOutputComponent("w","default","the real component of quaternion");
  keys.addOutputComponent("i","default","the i component of the quaternion");
  keys.addOutputComponent("j","default","the j component of the quaternion");
  keys.addOutputComponent("k","default","the k component of the quaternion");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
}

Quaternion::Quaternion(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  value(4),
  derivs(4),
  virial(4)
{
  for(unsigned i=0; i<4; ++i) derivs[i].resize(3);
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  if(atoms.size()!=3) error("Number of specified atoms should be 3");
  // please note, I do NO checking if these atoms are in the same molecule at all, so be careful in your input 

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  unsigned mode = getModeAndSetupValues( this );
  requestAtoms(atoms);
}

void Quaternion::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOMS",num,t);
  if( t.size()>0 ) aa->log.printf("  involving atoms %d %d %d\n",t[0].serial(),t[1].serial(),t[0].serial());
}

unsigned Quaternion::getModeAndSetupValues( ActionWithValue* av ) {
  // This sets up values that we can pass around in PLUMED
  av->addComponentWithDerivatives("w"); av->componentIsNotPeriodic("w");
  av->addComponentWithDerivatives("i"); av->componentIsNotPeriodic("i");
  av->addComponentWithDerivatives("j"); av->componentIsNotPeriodic("j");
  av->addComponentWithDerivatives("k"); av->componentIsNotPeriodic("k");
  return 0;
}

void Quaternion::calculate() {
  if(pbc) makeWhole();

  calculateCV( 0, masses, charges, getPositions(), value, derivs, virial, this );
  for(unsigned j=0;j<4;++j) {
      Value* valuej=getPntrToComponent(j);
      for(unsigned i=0;i<3;++i) setAtomsDerivatives(valuej,i,derivs[j][i] );
      setBoxDerivatives(valuej,virial[j]);
      valuej->set(value[j]);
  }
}

// calculator
void Quaternion::calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges,
                              const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs,
                              std::vector<Tensor>& virial, const ActionAtomistic* aa ) {
  //declarations
  Vector vec1_comp = delta( pos[0], pos[1] ); //components between atom 1 and 2
  Vector vec2_comp = delta( pos[0], pos[2] ); //components between atom 1 and 3 

////////////////////////////////////
////////x-vector calculations///////
  double magx = vec1_comp.modulo(); Vector xt = vec1_comp / magx;
  std::vector<Tensor> dx(3); double magx3= magx*magx*magx;
  dx[0](0,0) = ( -(vec1_comp[1]*vec1_comp[1]+vec1_comp[2]*vec1_comp[2])/magx3 );   // dx/dx
  dx[0](0,1) = (  vec1_comp[0]*vec1_comp[1]/magx3 );                 // dx/dy
  dx[0](0,2) = (  vec1_comp[0]*vec1_comp[2]/magx3 );                 // dx/dz
  dx[0](1,0) = (  vec1_comp[1]*vec1_comp[0]/magx3 );                 // dy/dx
  dx[0](1,1) = ( -(vec1_comp[0]*vec1_comp[0]+vec1_comp[2]*vec1_comp[2])/magx3 );   // dy/dy
  dx[0](1,2) = (  vec1_comp[1]*vec1_comp[2]/magx3 );
  dx[0](2,0) = (  vec1_comp[2]*vec1_comp[0]/magx3 );
  dx[0](2,1) = (  vec1_comp[2]*vec1_comp[1]/magx3 );
  dx[0](2,2) = ( -(vec1_comp[1]*vec1_comp[1]+vec1_comp[0]*vec1_comp[0])/magx3 );

  dx[1](0,0) = ( (vec1_comp[1]*vec1_comp[1]+vec1_comp[2]*vec1_comp[2])/magx3 );
  dx[1](0,1) = ( -vec1_comp[0]*vec1_comp[1]/magx3 );
  dx[1](0,2) = ( -vec1_comp[0]*vec1_comp[2]/magx3 );
  dx[1](1,0) = ( -vec1_comp[1]*vec1_comp[0]/magx3 );
  dx[1](1,1) = ( (vec1_comp[0]*vec1_comp[0]+vec1_comp[2]*vec1_comp[2])/magx3 );
  dx[1](1,2) = ( -vec1_comp[1]*vec1_comp[2]/magx3 );
  dx[1](2,0) = ( -vec1_comp[2]*vec1_comp[0]/magx3 );
  dx[1](2,1) = ( -vec1_comp[2]*vec1_comp[1]/magx3 );
  dx[1](2,2) = ( (vec1_comp[1]*vec1_comp[1]+vec1_comp[0]*vec1_comp[0])/magx3 ); 
  dx[2].zero();

/////////////////////////////////////
////////y-vector calculations////////
  //project vec2_comp on to vec1_comp
  //first do dot product of unormalized x and unormed y, divided by magnitude of x^2
  double dp = dotProduct( vec1_comp, vec2_comp ); double magx2=magx*magx; 
  std::vector<Vector> fac_derivs(3); double magx4=magx2*magx2, fac = dp/magx2; //fac meaning factor on front
  fac_derivs[0] = (-vec2_comp - vec1_comp)/magx2 + 2*dp*vec1_comp / magx4;
  fac_derivs[1] = (vec2_comp)/(magx2) - 2*dp*vec1_comp / magx4;
  fac_derivs[2] = (vec1_comp)/(magx2);
  //now multiply fac by unormed x, and subtract it from unormed y, then normalize
  Vector yt = vec2_comp - fac*vec1_comp; std::vector<Tensor> dy(3);
  dy[0](0,0) = -1 - fac_derivs[0][0]*vec1_comp[0] + fac;   // dx/dx
  dy[0](0,1) = -fac_derivs[0][0]*vec1_comp[1];             // dx/dy
  dy[0](0,2) = -fac_derivs[0][0]*vec1_comp[2];
  dy[0](1,0) = -fac_derivs[0][1]*vec1_comp[0];
  dy[0](1,1) = -1 - fac_derivs[0][1]*vec1_comp[1] + fac;
  dy[0](1,2) = -fac_derivs[0][1]*vec1_comp[2];
  dy[0](2,0) = -fac_derivs[0][2]*vec1_comp[0];
  dy[0](2,1) = -fac_derivs[0][2]*vec1_comp[1];
  dy[0](2,2) = -1 - fac_derivs[0][2]*vec1_comp[2] + fac;

  dy[1](0,0) = -fac_derivs[1][0]*vec1_comp[0] - fac;
  dy[1](0,1) = -fac_derivs[1][0]*vec1_comp[1];
  dy[1](0,2) = -fac_derivs[1][0]*vec1_comp[2];
  dy[1](1,0) = -fac_derivs[1][1]*vec1_comp[0];
  dy[1](1,1) = -fac_derivs[1][1]*vec1_comp[1] - fac;
  dy[1](1,2) = -fac_derivs[1][1]*vec1_comp[2];
  dy[1](2,0) = -fac_derivs[1][2]*vec1_comp[0];
  dy[1](2,1) = -fac_derivs[1][2]*vec1_comp[1];
  dy[1](2,2) = -fac_derivs[1][2]*vec1_comp[2] - fac;

  dy[2](0,0) = 1 - fac_derivs[2][0]*vec1_comp[0];
  dy[2](0,1) = -fac_derivs[2][0]*vec1_comp[1];
  dy[2](0,2) = -fac_derivs[2][0]*vec1_comp[2];
  dy[2](1,0) = -fac_derivs[2][1]*vec1_comp[0];
  dy[2](1,1) = 1 - fac_derivs[2][1]*vec1_comp[1];
  dy[2](1,2) = -fac_derivs[2][1]*vec1_comp[2];
  dy[2](2,0) = -fac_derivs[2][2]*vec1_comp[0];
  dy[2](2,1) = -fac_derivs[2][2]*vec1_comp[1];
  dy[2](2,2) = 1 - fac_derivs[2][2]*vec1_comp[2];
  //now normalize, and we have our y vector
  double magy = yt.modulo(); double imagy = 1/magy, magy3 = magy*magy*magy;
  Tensor multi_y; for(unsigned i=0; i<3; ++i) multi_y.setRow(i, yt[i]*yt); 
  yt = yt / magy; for(unsigned i=0; i<3; ++i) dy[i] = dy[i] / magy - matmul( dy[i], multi_y ) / magy3; 
/////////////////////////////////////
///////z-vector calculations/////////

  //simply cross them to make z
  Vector zt = crossProduct(xt,yt); std::vector<Tensor> dz(3);
  dz[0].setRow( 0, crossProduct( dx[0].getRow(0), yt ) + crossProduct( xt, dy[0].getRow(0) ) );
  dz[0].setRow( 1, crossProduct( dx[0].getRow(1), yt ) + crossProduct( xt, dy[0].getRow(1) ) );
  dz[0].setRow( 2, crossProduct( dx[0].getRow(2), yt ) + crossProduct( xt, dy[0].getRow(2) ) );

  dz[1].setRow( 0, crossProduct( dx[1].getRow(0), yt ) + crossProduct( xt, dy[1].getRow(0) ) );
  dz[1].setRow( 1, crossProduct( dx[1].getRow(1), yt ) + crossProduct( xt, dy[1].getRow(1) ) );
  dz[1].setRow( 2, crossProduct( dx[1].getRow(2), yt ) + crossProduct( xt, dy[1].getRow(2) ) );

  dz[2].setRow( 0, crossProduct( xt, dy[2].getRow(0) ) );
  dz[2].setRow( 1, crossProduct( xt, dy[2].getRow(1) ) );
  dz[2].setRow( 2, crossProduct( xt, dy[2].getRow(2) ) );

//the above 9 components form an orthonormal basis, centered on the molecule in question
//the rotation matrix is generally the inverse of this matrix, and in this case since it is 1) orthogonal and 2) its determinant is 1
//the inverse is simply the transpose
//we can now convert it to a quaternion, which will then be used to rotate the vector separating 2 molecules - > q1_bar * r_21 * q1, where r21 is the vector separating the atoms
//r21 is a quaternion with a 0 real part ^^^
//essentially this rotates the vector into the frame of molecule 1

//QUICKFIX
//double a,b,c;
//a=x[0];b=x[1];c=x[2];
//double d,e,f;
//d=y[0];e=y[1];f=y[2];
//double g,h,n;
//g=z[0];h=z[1];n=z[2];
//
//x[0]=a;x[1]=d;x[2]=g;
//y[0]=b;y[1]=e;y[2]=h;
//z[0]=c;z[1]=f;z[2]=n;
//Vector x(xt[0],xt[1],xt[2]);
//Vector y(yt[0],yt[1],yt[2]);
//Vector z(zt[0],zt[1],zt[2]);
//JAKE QUICKFIX
  Vector x(xt[0],yt[0],zt[0]);
  Vector y(xt[1],yt[1],zt[1]);
  Vector z(xt[2],yt[2],zt[2]);

  double tr = x[0] + y[1] + z[2]; //trace of the rotation matrix

//aa->log.printf("%f %f %f\n%f %f %f\n%f %f %f\n",x[0],x[1],x[2],y[0],y[1],y[2],z[0],z[1],z[2]);
//this converts to quaternion
  std::vector<Vector> dS(3); 
  if (tr > 0.000001) { 
    double S = sqrt(tr+1.0) * 2; // S=4*qw
    for(unsigned i=0; i<3; ++i) dS[i] = (2/S)*(dx[i].getCol(0) + dy[i].getCol(1) + dz[i].getCol(2));  
    vals[0] = 0.25 * S; for(unsigned i=0; i<3; ++i) derivs[0][i] =0.25*dS[i];
    vals[1] = (z[1] - y[2]) / S;
    for(unsigned i=0; i<3; ++i) derivs[1][i] = (1/S)*(dy[i].getCol(2) - dz[i].getCol(1)) - (vals[1]/S)*dS[i];
    vals[2] = (x[2] - z[0]) / S; 
    for(unsigned i=0; i<3; ++i) derivs[2][i] = (1/S)*(dz[i].getCol(0) - dx[i].getCol(2)) - (vals[2]/S)*dS[i]; 
    vals[3] = (y[0] - x[1]) / S; 
    for(unsigned i=0; i<3; ++i) derivs[3][i] = (1/S)*(dx[i].getCol(1) - dy[i].getCol(0)) - (vals[3]/S)*dS[i];
  } else if ((x[0] > y[1])&(x[0] > z[2])) { 
    float S = sqrt(1.0 + x[0] - y[1] - z[2]) * 2; // S=4*qx 
    for(unsigned i=0; i<3; ++i) dS[i] = (2/S)*(dx[i].getCol(0) - dy[i].getCol(1) - dz[i].getCol(2));
    vals[0] = (z[1] - y[2]) / S;
    for(unsigned i=0; i<3; ++i) derivs[0][i] = (1/S)*(dz[i].getCol(1) - dy[i].getCol(2)) - (vals[0]/S)*dS[i];
    vals[1] = 0.25 * S; for(unsigned i=0; i<3; ++i) derivs[1][i] =0.25*dS[i];
    vals[2] = (x[1] + y[0]) / S; 
    for(unsigned i=0; i<3; ++i) derivs[2][i] = (1/S)*(dx[i].getCol(1) + dy[i].getCol(0)) - (vals[2]/S)*dS[i];
    vals[3] = (x[2] + z[0]) / S; 
    for(unsigned i=0; i<3; ++i) derivs[3][i] = (1/S)*(dx[i].getCol(2) + dz[i].getCol(0)) - (vals[3]/S)*dS[i];
  } else if (y[1] > z[2]) { 
    float S = sqrt(1.0 + y[1] - x[0] - z[2]) * 2; // S=4*qy
    for(unsigned i=0; i<3; ++i) dS[i] = (2/S)*( -dx[i].getCol(0) + dy[i].getCol(1) - dz[i].getCol(2)); 
    vals[0] = (x[2] - z[0]) / S;
    for(unsigned i=0; i<3; ++i) derivs[0][i] = (1/S)*(dx[i].getCol(2) - dz[i].getCol(0)) - (vals[0]/S)*dS[i]; 
    vals[1] = (x[1] + y[0]) / S; 
    for(unsigned i=0; i<3; ++i) derivs[1][i] = (1/S)*(dx[i].getCol(1) + dy[i].getCol(0)) - (vals[1]/S)*dS[i];
    vals[2] = 0.25 * S; for(unsigned i=0; i<3; ++i) derivs[2][i] =0.25*dS[i];
    vals[3] = (y[2] + z[1]) / S; 
    for(unsigned i=0; i<3; ++i) derivs[3][i] = (1/S)*(dy[i].getCol(2) + dz[i].getCol(1)) - (vals[3]/S)*dS[i]; 
  } else {
    float S = sqrt(1.0 + z[2] - x[0] - y[1]) * 2; // S=4*qz
    for(unsigned i=0; i<3; ++i) dS[i] = (2/S)*(-dx[i].getCol(0) - dy[i].getCol(1) + dz[i].getCol(2)); 
    vals[0] = (y[0] - x[1]) / S;
    for(unsigned i=0; i<3; ++i) derivs[0][i] = (1/S)*(dy[i].getCol(0) - dx[i].getCol(1)) - (vals[0]/S)*dS[i];
    vals[1] = (x[2] + z[0]) / S;
    for(unsigned i=0; i<3; ++i) derivs[1][i] = (1/S)*(dx[i].getCol(2) + dz[i].getCol(0)) - (vals[1]/S)*dS[i]; 
    vals[2] = (y[2] + z[1]) / S;
    for(unsigned i=0; i<3; ++i) derivs[2][i] = (1/S)*(dy[i].getCol(2) + dz[i].getCol(1)) - (vals[2]/S)*dS[i]; 
    vals[3] = 0.25 * S; for(unsigned i=0; i<3; ++i) derivs[3][i] =0.25*dS[i];
  }
  //JAKE QUICKFIX
  //if (vals[0] < 0) {vals[0]*=-1;vals[1]*=-1;vals[2]*=-1;vals[3]*=-1;}
  setBoxDerivativesNoPbc( pos, derivs, virial );
//worked but derivatvies are fucked

//debugging vector calcs
//  aa->log.printf("%f ", derivs[0]);
//  aa->log.printf("%f ", derivs[1]);
//  aa->log.printf("%f ", derivs[2]);
//  aa->log.printf("%f\n", derivs[3]);
//log.printf("%f ", y2);
//log.printf("%f\n", y3);
//log.printf("%f ",z1);
//log.printf("%f ",z2);
//log.printf("%f\n",z3);

}

}
}



