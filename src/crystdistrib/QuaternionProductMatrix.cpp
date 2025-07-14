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

Calculate the outer product matrix from two vectors of quaternions. Quaternions are made of four numbers, one real number, and three imaginary numbers \textit{i}, \textit{j}, and \textit{k}. The imaginary numbers are not commutative:

\[ij =k\] \[jk=i\] \[ki=j\] \[ji = -k\] \[ik = -j\] \[kj = -i\] \[ii = jj = kk = -1\]

Thus multiplying two quaternions $\mathbf{q_1} = w_1 + x_1i + y_1j + z_1k  $ and $\mathbf{q_2} = w_2 + x_2i + y_2j + z_2k$ yields the following four components:

\[w_{12} = w_1w_2 - x_1x_2 - y_1y_2 -z_1z_2\]
\[x_{12}i = (w_1x_2 + x_1w_2 + y_1z_2 - z_1y_2  )i\]
\[y_{12}j = (w_1y_2 - x_1z_2 + y_1w_2 + z_1x_2)j\]
\[z_{12}k = (w_1z_2 + x_1y_2 - y_1x_2 + z_1w_2)k\]

Quaternions can also be multiplied by a real number, and the conjugate $\mathbf{\overline{q}} = w -xi - yj - zk$ is analogous to complex numbers.


This action takes two equal sized vectors of quaternions of length $n$, and returns four $n\times n$ outer product matrices, corresponding to components w, x, y, and z in the above example. These matrices are real numbers, and will not behave with the usual nuances of quaternion algebra in any CUSTOMS, or other actions, that will need to be accounted for manually. The vectors of quaternions can be the same vectors, or different. Note, since the QUATERNION action returns unit quaternions, all quaternions returned here will also be unit quaternions.

```plumed
#in a system of 4 molecules, each with 3 atoms
#define quaternions 1-4
quats12: QUATERNION ATOMS1=1,2,3 ATOMS2=4,5,6
quats34: QUATERNION ATOMS1=7,8,9 ATOMS2=10,11,12

#take the outer product of 1,2 and 3,4
qp: QUATERNION_PRODUCT_MATRIX ARG=quats12.*,quats34.* #notice the asterisk
#the components need to be passed in the order w1,x1,y1,z1,w2,x2,y2,z2
#now use in an order parameter
bpw: CUSTOM ARG=qp.* VAR=w,x,y,z FUNC=exp(w/2+x/2+y/2+z/2) PERIODIC=NO
```

A common strategy when using this action is to multiply a function of elements of the quaternion product matrix by a contact matrix as shown below:

```plumed
# Calculate some quaternions
quats12: QUATERNION ATOMS1=1,2,3 ATOMS2=4,5,6 ATOMS3=7,8,9 ATOMS4=10,11,12
# Calculate a contact matrix
cmap: CONTACT_MATRIX GROUP=1,4,7,10 SWITCH={RATIONAL R_0=0.2 D_MAX=0.7}
# Calculate the quaternion product matrix
qp: QUATERNION_PRODUCT_MATRIX ARG=quats12.*,quats12.* MASK=cmap
# Now calculate a function of the quaternion elements (normally something more complicated than the following is performed
func: CUSTOM ARG=cmap,qp.* VAR=x,w,i,j,k FUNC=x*(w+i+j+k) PERIODIC=NO
# Now sum the rows of the above matrix
ones: ONES SIZE=4
op: MATRIX_VECTOR_PRODUCT ARG=func,ones
# And output the values of the order parameter
DUMPATOMS ATOMS=1,4,7,10 ARG=op FILE=func.yxz
```

Notice how the MASK keyword is used in the CUSTOM in the input to QUATERNION_PRODUCT_MATRIX to prevent PLUMED from calculating the
elements of the QUATERNION_PRODUCT_MATRIX that will be multiplied by the elements of that matrix `cmap` that is calculated in the
CONTACT_MATRIX action when the calculations in the CUSTOM action with label `func` are performed.  This is the same procedure that is performed
in [TORSIONS_MATRIX](TORSIONS_MATRIX.md) and [MATRIX_PRODUCT](MATRIX_PRODUCT.md) to avoid computational expense.  This procedure is automatically
performed when you use the [ROPS](ROPS.md) shortcut.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystdistrib {

class QuaternionProductMatrix {
public:
  static void registerKeywords( Keywords& keys );
  void setup( const std::vector<std::size_t>& shape, const std::string& func, matrixtools::OuterProductBase<QuaternionProductMatrix>* action );
  static void calculate( bool noderiv, const QuaternionProductMatrix& actdata, View<double> vals, MatrixElementOutput& output );
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

void QuaternionProductMatrix::calculate( bool noderiv, const QuaternionProductMatrix& actdata, View<double> vals, MatrixElementOutput& output ) {
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
