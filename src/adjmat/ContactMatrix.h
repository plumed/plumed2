/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#ifndef __PLUMED_adjmat_ContactMatrix_h
#define __PLUMED_adjmat_ContactMatrix_h

#include "AdjacencyMatrixBase.h"
#include "tools/SwitchingFunction.h"
namespace PLMD {
namespace adjmat {

struct ContactMatrix {
#ifdef __PLUMED_HAS_OPENACC
  SwitchingFunctionAccelerable switchingFunction;
#else
  SwitchingFunction switchingFunction;
#endif
  static void registerKeywords( Keywords& keys );
  template <typename myPTM>
  void parseInput( AdjacencyMatrixBase<ContactMatrix,myPTM>* action );
  static void calculateWeight( const ContactMatrix& data,
                               const AdjacencyMatrixInput& input,
                               MatrixOutput output );
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(this[0:1])
  }
#endif //__PLUMED_HAS_OPENACC
};

void ContactMatrix::registerKeywords( Keywords& keys ) {
  keys.setDisplayName("CONTACT_MATRIX");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
}

template <typename myPTM>
void ContactMatrix::parseInput( AdjacencyMatrixBase<ContactMatrix,myPTM>* action ) {
  std::string errors;
  std::string swinput;
  action->parse("SWITCH",swinput);
  if( swinput.length()>0 ) {
    switchingFunction.set( swinput, errors );
    if( errors.length()!=0 ) {
      action->error("problem reading switching function description " + errors);
    }
  } else {
    int nn=0;
    int mm=0;
    double r_0=-1.0;
    double d_0= 0.0;
    action->parse("NN",nn);
    action->parse("MM",mm);
    action->parse("R_0",r_0);
    action->parse("D_0",d_0);
    if( r_0<0.0 ) {
      action->error("you must set a value for R_0");
    }
    switchingFunction.set(nn,mm,r_0,d_0);
  }
  // And set the link cell cutoff
  action->log.printf("  switching function cutoff is %s \n",switchingFunction.description().c_str() );
  action->setLinkCellCutoff( true, switchingFunction.get_dmax() );
}

void ContactMatrix::calculateWeight( const ContactMatrix& data,
                                     const AdjacencyMatrixInput& input,
                                     MatrixOutput output ) {
  const double mod2 = input.pos.modulo2();
  if( mod2<epsilon ) {
    return;  // Atoms can't be bonded to themselves
  }
  double dfunc;
  output.val[0] = data.switchingFunction.calculateSqr( mod2, dfunc );
  if( output.val[0]<epsilon ) {
    output.val[0] = 0.0;
    return;
  }
  if( input.noderiv ) {
    return;
  }
  const Vector v { (-dfunc)*input.pos[0],
                   (-dfunc)*input.pos[1],
                   (-dfunc)*input.pos[2] };
  output.deriv[0] = v[0];
  output.deriv[1] = v[1];
  output.deriv[2] = v[2];
  output.deriv[3] =-v[0];
  output.deriv[4] =-v[1];
  output.deriv[5] =-v[2];

  output.assignOuterProduct( 6, v, input.pos);
}
}
}
#endif
