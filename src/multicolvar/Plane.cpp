/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

class Plane : public MultiColvarBase {
private:
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void createVectorNormInput( const std::string& ilab, const std::string& olab, 
                                     const std::string& vlab, std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit Plane(const ActionOptions&);
// active methods:
  void compute( const unsigned& index, AtomValuePack& myatoms ) const ; 
};

PLUMED_REGISTER_ACTION(Plane,"PLANES")
PLUMED_REGISTER_SHORTCUT(Plane,"PLANES")

void Plane::shortcutKeywords( Keywords& keys ) {
  keys.addFlag("VMEAN",false,"calculate the norm of the mean vector.");
  keys.addOutputComponent("_vmean","VMEAN","the norm of the mean vector");
  keys.addFlag("VSUM",false,"calculate the norm of the sum of all the vectors");
  keys.addOutputComponent("_vsum","VSUM","the norm of the mean vector");
} 

void Plane::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                            const std::map<std::string,std::string>& keys,
                            std::vector<std::vector<std::string> >& actions ) {
  std::vector<std::string> pl_input; pl_input.push_back( lab + ":" );
  pl_input.push_back("PLANES"); for(unsigned i=1;i<words.size();++i) pl_input.push_back( words[i] );
  actions.push_back( pl_input );
  if( keys.count("VMEAN") ) {
      std::vector<std::string> x_input; x_input.push_back(lab + "_xs:"); 
      x_input.push_back("COMBINE"); x_input.push_back("ARG=" + lab + ".x");
      x_input.push_back("NORMALIZE"); x_input.push_back("PERIODIC=NO");
      actions.push_back(x_input);
      std::vector<std::string> y_input; y_input.push_back(lab + "_ys:"); 
      y_input.push_back("COMBINE"); y_input.push_back("ARG=" + lab + ".y");
      y_input.push_back("NORMALIZE"); y_input.push_back("PERIODIC=NO");
      actions.push_back(y_input);
      std::vector<std::string> z_input; z_input.push_back(lab + "_zs:"); 
      z_input.push_back("COMBINE"); z_input.push_back("ARG=" + lab + ".z");
      z_input.push_back("NORMALIZE"); z_input.push_back("PERIODIC=NO");
      actions.push_back(z_input);
      // Now calculate the total length of the vector
      createVectorNormInput( lab, lab + "_vmean", "s", actions );
  }  
  if( keys.count("VSUM") ) {
      std::vector<std::string> x_input; x_input.push_back(lab + "_xz:"); 
      x_input.push_back("COMBINE"); x_input.push_back("ARG=" + lab + ".x");
      x_input.push_back("PERIODIC=NO"); actions.push_back(x_input);
      std::vector<std::string> y_input; y_input.push_back(lab + "_yz:");
      y_input.push_back("COMBINE"); y_input.push_back("ARG=" + lab + ".y");
      y_input.push_back("PERIODIC=NO"); actions.push_back(y_input);
      std::vector<std::string> z_input; z_input.push_back(lab + "_zz:");
      z_input.push_back("COMBINE"); z_input.push_back("ARG=" + lab + ".z");
      z_input.push_back("PERIODIC=NO"); actions.push_back(z_input);
      // Now calculate the total length of the vector
      createVectorNormInput( lab, lab + "_vsum", "z", actions );
  }
}

void Plane::createVectorNormInput( const std::string& ilab, const std::string& olab, 
                                   const std::string& vlab, std::vector<std::vector<std::string> >& actions ) {
  std::vector<std::string> norm_input; norm_input.push_back(olab + "2:");
  norm_input.push_back("COMBINE"); std::string powstr="POWERS=2";
  norm_input.push_back("ARG1=" + ilab + "_x" + vlab );
  norm_input.push_back("ARG2=" + ilab + "_y" + vlab );
  norm_input.push_back("ARG3=" + ilab + "_z" + vlab );
  norm_input.push_back("POWERS=2,2,2"); norm_input.push_back("PERIODIC=NO"); actions.push_back( norm_input );
  std::vector<std::string> sqrt_input; sqrt_input.push_back(olab + ":"); sqrt_input.push_back("MATHEVAL");
  sqrt_input.push_back("ARG1=" + olab + "2"); sqrt_input.push_back( "FUNC=sqrt(x)" );
  sqrt_input.push_back("PERIODIC=NO"); actions.push_back( sqrt_input );
}


void Plane::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.addOutputComponent("x","default","the x-component of the vector that is normal to the plane containing the atoms");
  keys.addOutputComponent("y","default","the y-component of the vector that is normal to the plane containing the atoms");
  keys.addOutputComponent("z","default","the z-component of the vector that is normal to the plane containing the atoms");
}

Plane::Plane(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  if(getNumberOfAtomsInEachCV()!=3 && getNumberOfAtomsInEachCV()!=4){
     error("Number of specified atoms should be 3 or 4");
  }
  checkRead();
  addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
  addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
}


// calculator
void Plane::compute( const unsigned& index, AtomValuePack& myatoms ) const {

  Vector d1, d2, cp;
  if( myatoms.getNumberOfAtoms()==3 ) {
    d1=delta( myatoms.getPosition(1), myatoms.getPosition(0) );
    d2=delta( myatoms.getPosition(1), myatoms.getPosition(2) );
  } else {
    d1=delta( myatoms.getPosition(1), myatoms.getPosition(0) );
    d2=delta( myatoms.getPosition(2), myatoms.getPosition(3) );
  }
  cp = crossProduct( d1, d2 );

  myatoms.addAtomsDerivatives( 0, 0, crossProduct( Vector(-1.0,0,0), d2 ));
  if( myatoms.getNumberOfAtoms()==3 ) {
    myatoms.addAtomsDerivatives( 0, 1, crossProduct( Vector(+1.0,0,0), d2 ) + crossProduct( Vector(-1.0,0,0), d1 ));
    myatoms.addAtomsDerivatives( 0, 2, crossProduct( Vector(+1.0,0,0), d1 ));
  } else {
    myatoms.addAtomsDerivatives( 0, 1, crossProduct( Vector(+1.0,0,0), d2 ));
    myatoms.addAtomsDerivatives( 0, 2, crossProduct( Vector(-1.0,0,0), d1 ));
    myatoms.addAtomsDerivatives( 0, 3, crossProduct( Vector(+1.0,0,0), d1 ));
  }
  myatoms.addBoxDerivatives( 0, Tensor(d1,crossProduct(Vector(+1.0,0,0), d2)) + Tensor( d2, crossProduct(Vector(-1.0,0,0), d1)) );
  myatoms.setValue( 0, cp[0] );

  myatoms.addAtomsDerivatives( 1, 0, crossProduct( Vector(0,-1.0,0), d2 ));
  if( myatoms.getNumberOfAtoms()==3 ) {
    myatoms.addAtomsDerivatives( 1, 1, crossProduct( Vector(0,+1.0,0), d2 ) + crossProduct( Vector(0,-1.0,0), d1 ));
    myatoms.addAtomsDerivatives( 1, 2, crossProduct( Vector(0,+1.0,0), d1 ));
  } else {
    myatoms.addAtomsDerivatives( 1, 1, crossProduct( Vector(0,+1.0,0), d2 ));
    myatoms.addAtomsDerivatives( 1, 2, crossProduct( Vector(0,-1.0,0), d1 ));
    myatoms.addAtomsDerivatives( 1, 3, crossProduct( Vector(0,+1.0,0), d1 ));
  }
  myatoms.addBoxDerivatives( 1, Tensor(d1,crossProduct(Vector(0,+1.0,0), d2)) + Tensor( d2, crossProduct(Vector(0,-1.0,0), d1)) );
  myatoms.setValue( 1, cp[1] );

  myatoms.addAtomsDerivatives( 2, 0, crossProduct( Vector(0,0,-1.0), d2 ));
  if( myatoms.getNumberOfAtoms()==3 ) {
    myatoms.addAtomsDerivatives( 2, 1, crossProduct( Vector(0,0,+1.0), d2 ) + crossProduct( Vector(0,0,-1.0), d1 ));
    myatoms.addAtomsDerivatives( 2, 2, crossProduct( Vector(0,0,+1.0), d1 ));
  } else {
    myatoms.addAtomsDerivatives( 2, 1, crossProduct( Vector(0,0,-1.0), d2 ));
    myatoms.addAtomsDerivatives( 2, 2, crossProduct( Vector(0,0,-1.0), d1 ));
    myatoms.addAtomsDerivatives( 2, 3, crossProduct( Vector(0,0,+1.0), d1 ));
  }
  myatoms.addBoxDerivatives( 2, Tensor(d1,crossProduct(Vector(0,0,+1.0), d2)) + Tensor( d2, crossProduct(Vector(0,0,-1.0), d1)) );
  myatoms.addValue( 2, cp[2] );

}

}
}



