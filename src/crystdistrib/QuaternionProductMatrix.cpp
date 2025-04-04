/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "matrixtools/OuterProduct.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR QUATERNION_PRODUCT_MATRIX
/*
Calculate the outer product matrix from two vectors of quaternions

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystdistrib {

class QuaternionProductMatrix {
public:
  static void registerKeywords( Keywords& keys );
  void setup( const std::vector<std::size_t>& shape, const std::string& func, matrixtools::OuterProductBase<QuaternionProductMatrix>* action );
  static void calculate( bool noderiv, const QuaternionProductMatrix& actdata, const std::vector<double>& vals, MatrixElementOutput& output );
};

typedef matrixtools::OuterProductBase<QuaternionProductMatrix> qpop;
PLUMED_REGISTER_ACTION(qpop,"QUATERNION_PRODUCT_MATRIX")

void QuaternionProductMatrix::registerKeywords( Keywords& keys ) {
  keys.addInputKeyword("compulsory","ARG","vector","the labels of the quaternion vectors that you are outer product of");
  keys.addInputKeyword("optional","MASK","matrix","a matrix that is used to used to determine which elements of the output matrix to compute");
  keys.addOutputComponent("w","default","matrix","the real component of quaternion");
  keys.addOutputComponent("i","default","matrix","the i component of the quaternion");
  keys.addOutputComponent("j","default","matrix","the j component of the quaternion");
  keys.addOutputComponent("k","default","matrix","the k component of the quaternion");
}

void QuaternionProductMatrix::setup( const std::vector<std::size_t>& shape, const std::string& func, matrixtools::OuterProductBase<QuaternionProductMatrix>* action ) {
  unsigned nargs=action->getNumberOfArguments();
  if( action->getNumberOfMasks()>0 ) {
    nargs = nargs - action->getNumberOfMasks();
  }
  if( nargs!=8 ) {
    action->error("should be eight arguments to this action.  Four quaternions for each set of atoms.  You can repeat actions");
  }
  for(unsigned i=0; i<8; ++i) {
    Value* myarg=action->getPntrToArgument(i);
    if( myarg->getRank()!=1 ) {
      action->error("all arguments to this action should be vectors");
    }
    if( (myarg->getPntrToAction())->getName()!="QUATERNION_VECTOR" ) {
      action->error("all arguments to this action should be quaternions");
    }
    std::string mylab=(action->getPntrToArgument(i))->getName();
    std::size_t dot=mylab.find_first_of(".");
    if( (i==0 || i==4) && mylab.substr(dot+1)!="w" ) {
      action->error("quaternion arguments are in wrong order");
    }
    if( (i==1 || i==5) && mylab.substr(dot+1)!="i" ) {
      action->error("quaternion arguments are in wrong order");
    }
    if( (i==2 || i==6) && mylab.substr(dot+1)!="j" ) {
      action->error("quaternion arguments are in wrong order");
    }
    if( (i==3 || i==7) && mylab.substr(dot+1)!="k" ) {
      action->error("quaternion arguments are in wrong order");
    }
  }
  action->addComponent( "w", shape );
  action->componentIsNotPeriodic("w");
  action->addComponent( "i", shape );
  action->componentIsNotPeriodic("i");
  action->addComponent( "j", shape );
  action->componentIsNotPeriodic("j");
  action->addComponent( "k", shape );
  action->componentIsNotPeriodic("k");
}

void QuaternionProductMatrix::calculate( bool noderiv, const QuaternionProductMatrix& actdata, const std::vector<double>& vals, MatrixElementOutput& output ) {
  std::vector<double> quat1(4), quat2(4);

  for(unsigned i=0; i<4; ++i) {
    // Retrieve the first quaternion
    quat1[i] = vals[i];
    // Retrieve the second quaternion
    quat2[i] = vals[4+i];
  }

  //make q1 the conjugate
  quat1[1] *= -1;
  quat1[2] *= -1;
  quat1[3] *= -1;

  double pref=1;
  double pref2=1;
  double conj=1;
  //real part of q1*q2
  for(unsigned i=0; i<4; ++i) {
    if( i>0 ) {
      pref=-1;
      pref2=-1;
    }
    output.values[0] += pref*quat1[i]*quat2[i];
    if( noderiv ) {
      continue ;
    }
    if (i>0) {
      conj=-1;
    }
    output.derivs[0][i] = conj*pref*quat2[i];
    output.derivs[0][4+i] = pref2*quat1[i];
  }
  //i component
  pref=1;
  conj=1;
  pref2=1;
  for(unsigned i=0; i<4; i++) {
    if(i==3) {
      pref=-1;
    } else {
      pref=1;
    }
    if(i==2) {
      pref2=-1;
    } else {
      pref2=1;
    }
    output.values[1] += pref*quat1[i]*quat2[(5-i)%4];
    if( noderiv ) {
      continue ;
    }
    if (i>0) {
      conj=-1;
    }
    output.derivs[1][i] = conj*pref*quat2[(5-i)%4];
    output.derivs[1][4+i] = pref2*quat1[(5-i)%4];
  }

  //j component
  pref=1;
  conj=1;
  pref2=1;
  for (unsigned i=0; i<4; i++) {
    if(i==1) {
      pref=-1;
    } else {
      pref=1;
    }
    if (i==3) {
      pref2=-1;
    } else {
      pref2=1;
    }
    output.values[2] += pref*quat1[i]*quat2[(i+2)%4];
    if( noderiv ) {
      continue ;
    }
    if (i>0) {
      conj=-1;
    }
    output.derivs[2][i] = conj*pref*quat2[(i+2)%4];
    output.derivs[2][4+i] = pref2*quat1[(i+2)%4];
  }

  //k component
  pref=1;
  conj=1;
  pref2=1;
  for (unsigned i=0; i<4; i++) {
    if(i==2) {
      pref=-1;
    } else {
      pref=1;
    }
    if(i==1) {
      pref2=-1;
    } else {
      pref2=1;
    }
    output.values[3] += pref*quat1[i]*quat2[(3-i)];
    if( noderiv ) {
      continue ;
    }
    if (i>0) {
      conj=-1;
    }
    output.derivs[3][i] = conj*pref*quat2[3-i];
    output.derivs[3][4+i] = pref2*quat1[3-i];
  }
}

}
}
