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
  for(unsigned i=0; i<4; ++i) derivs[i].resize(2);
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
      
  double x1, x2, x3; //vector 1
  double y1, y2, y3; //vector 2
  double z1, z2, z3; //vector 3

////////////////////////////////////
////////x-vector calculations///////
  double vec1_magnitude = sqrt((vec1_comp[0]*vec1_comp[0])+(vec1_comp[1]*vec1_comp[1])+(vec1_comp[2]*vec1_comp[2])); //normalization factor for atom 2 - atom 1

  x1 = vec1_comp[0]/vec1_magnitude;
  x2 = vec1_comp[1]/vec1_magnitude;
  x3 = vec1_comp[2]/vec1_magnitude;

/////////////////////////////////////
////////y-vector calculations////////
  //project vec2_comp on to vec1_comp
  //first do dot product of unormalized x and unormed y, divided by magnitude of x^2
  double fac = ((vec1_comp[0]*vec2_comp[0])+(vec1_comp[1]*vec2_comp[1])+(vec1_comp[2]*vec2_comp[2]))/(vec1_magnitude*vec1_magnitude); //fac meaning factor on front
  //now multiply fac by unormed x, and subtract it from unormed y, then normalize
  y1 = vec2_comp[0] - fac*vec1_comp[0];
  y2 = vec2_comp[1] - fac*vec1_comp[1];
  y3 = vec2_comp[2] - fac*vec1_comp[2];
  //now normalize, and we have our y vector
  double y_magnitude = sqrt((y1*y1)+(y2*y2)+(y3*y3));
  y1 = y1/y_magnitude;
  y2 = y2/y_magnitude;
  y3 = y3/y_magnitude;
/////////////////////////////////////
///////z-vector calculations/////////

  //simply cross them to make z
  Vector cross_xy = crossProduct(Vector(x1,x2,x3),Vector(y1,y2,y3));

  z1 = cross_xy[0];
  z2 = cross_xy[1];
  z3 = cross_xy[2];

//the above 9 components form an orthonormal basis, centered on the molecule in question
//the rotation matrix is generally the inverse of this matrix, and in this case since it is 1) orthogonal and 2) its determinant is 1
//the inverse is simply the transpose
//we can now convert it to a quaternion, which will then be used to rotate the vector separating 2 molecules - > q1_bar * r_21 * q1, where r21 is the vector separating the atoms
//r21 is a quaternion with a 0 real part ^^^
//essentially this rotates the vector into the frame of molecule 1
  double tr = x1 + y2 + z3; //trace of the rotation matrix

//this converts to quaternion
  if (tr > 0) { 
    double S = sqrt(tr+1.0) * 2; // S=4*qw 
    vals[0] = 0.25 * S;
    vals[1] = (z2 - y3) / S;
    vals[2] = (x3 - z1) / S; 
    vals[3] = (y1 - x2) / S; 
  } else if ((x1 > y2)&(x1 > z3)) { 
    float S = sqrt(1.0 + x1 - y2 - z3) * 2; // S=4*qx 
    vals[0] = (z2 - y3) / S;
    vals[1] = 0.25 * S;
    vals[2] = (x2 + y1) / S; 
    vals[3] = (x3 + z1) / S; 
  } else if (y2 > z3) { 
    float S = sqrt(1.0 + y2 - x1 - z3) * 2; // S=4*qy
    vals[0] = (x3 - z1) / S;
    vals[1] = (x2 + y1) / S; 
    vals[2] = 0.25 * S;
    vals[3] = (y3 + z2) / S; 
  } else { 
    float S = sqrt(1.0 + z3 - x1 - y2) * 2; // S=4*qz
    vals[0] = (y1 - x2) / S;
    vals[1] = (x3 + z1) / S;
    vals[2] = (y3 + z2) / S;
    vals[3] = 0.25 * S;
  }
//debugging vector calcs
//log.printf("%f ", x1);
//log.printf("%f ", x2);
//log.printf("%f\n", x3);
//log.printf("%f ", y1);
//log.printf("%f ", y2);
//log.printf("%f\n", y3);
//log.printf("%f ",z1);
//log.printf("%f ",z2);
//log.printf("%f\n",z3);

}

}
}



