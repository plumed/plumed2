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
#include "core/ActionWithVector.h"
#include "core/ActionRegister.h"
#include "tools/Torsion.h"


#include <iostream>

namespace PLMD {
namespace crystdistrib {

//+PLUMEDOC MCOLVAR QUATERNION_BOND_PRODUCT_MATRIX
/*
Calculate the product between a matrix of quaternions and the bonds

\par Examples

*/
//+ENDPLUMEDOC

class QuaternionBondProductMatrix : public ActionWithVector {
private:
  std::vector<unsigned> active_tasks;
//  const Vector4d& rightMultiply(Tensor4d&, Vector4d&);
public:
  static void registerKeywords( Keywords& keys );
  explicit QuaternionBondProductMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  void prepare() override ;
  void calculate() override ;
  std::vector<unsigned>& getListOfActiveTasks( ActionWithVector* action ) override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
};

PLUMED_REGISTER_ACTION(QuaternionBondProductMatrix,"QUATERNION_BOND_PRODUCT_MATRIX")


//const Vector4d& QuaternionBondMatrix::rightMultiply(Tensor4d& pref, Vector4d& quat) {
//  Vector4d temp;
//  int sumTemp;
//  for (int i=0; i<4; i++){ //rows
//    sumTemp=0;
//    for (int j=0; j<4; j++){ //cols
//      sumTemp+=pref(i,j)*quat[j];
//    }
//    temp[i]=sumTemp;
//  }
//  return temp;
//}




void QuaternionBondProductMatrix::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","vector/matrix","this action takes 8 arguments.  The first four should be the w,i,j and k components of a quaternion vector.  The second four should be contact matrix and the matrices should be the x, y and z components of the bond vectors");
  keys.addOutputComponent("w","default","matrix","the real component of quaternion");
  keys.addOutputComponent("i","default","matrix","the i component of the quaternion");
  keys.addOutputComponent("j","default","matrix","the j component of the quaternion");
  keys.addOutputComponent("k","default","matrix","the k component of the quaternion");
}

QuaternionBondProductMatrix::QuaternionBondProductMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao) {
  if( getNumberOfArguments()!=8 ) {
    error("should be eight arguments to this action, 4 quaternion components and 4 matrices");
  }
  unsigned nquat = getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=0; i<4; ++i) {
    Value* myarg=getPntrToArgument(i);
    if( myarg->getRank()!=1 ) {
      error("first four arguments to this action should be vectors");
    }
    if( (myarg->getPntrToAction())->getName()!="QUATERNION_VECTOR" ) {
      error("first four arguments to this action should be quaternions");
    }
    std::string mylab=getPntrToArgument(i)->getName();
    std::size_t dot=mylab.find_first_of(".");
    if( i==0 && mylab.substr(dot+1)!="w" ) {
      error("quaternion arguments are in wrong order");
    }
    if( i==1 && mylab.substr(dot+1)!="i" ) {
      error("quaternion arguments are in wrong order");
    }
    if( i==2 && mylab.substr(dot+1)!="j" ) {
      error("quaternion arguments are in wrong order");
    }
    if( i==3 && mylab.substr(dot+1)!="k" ) {
      error("quaternion arguments are in wrong order");
    }
  }
  std::vector<unsigned> shape( getPntrToArgument(4)->getShape() );
  for(unsigned i=4; i<8; ++i) {
    Value* myarg=getPntrToArgument(i);
    if( myarg->getRank()!=2 ) {
      error("second four arguments to this action should be matrices");
    }
    if( myarg->getShape()[0]!=shape[0] || myarg->getShape()[1]!=shape[1] ) {
      error("matrices should all have the same shape");
    }
    if( myarg->getShape()[0]!=nquat ) {
      error("number of rows in matrix should equal number of input quaternions");
    }
    std::string mylab=getPntrToArgument(i)->getName();
    std::size_t dot=mylab.find_first_of(".");
    if( i==5 && mylab.substr(dot+1)!="x" ) {
      error("quaternion arguments are in wrong order");
    }
    if( i==6 && mylab.substr(dot+1)!="y" ) {
      error("quaternion arguments are in wrong order");
    }
    if( i==7 && mylab.substr(dot+1)!="z" ) {
      error("quaternion arguments are in wrong order");
    }
  }
  addComponent( "w", shape );
  componentIsNotPeriodic("w");
  addComponent( "i", shape );
  componentIsNotPeriodic("i");
  addComponent( "j", shape );
  componentIsNotPeriodic("j");
  addComponent( "k", shape );
  componentIsNotPeriodic("k");
}

unsigned QuaternionBondProductMatrix::getNumberOfDerivatives() {
  unsigned nder=0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    nder += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  return nder;
}

void QuaternionBondProductMatrix::prepare() {
  ActionWithVector::prepare();
  active_tasks.resize(0);
}

std::vector<unsigned>& QuaternionBondProductMatrix::getListOfActiveTasks( ActionWithVector* action ) {
  if( active_tasks.size()>0 ) {
    return active_tasks;
  }

  Value* myarg = getPntrToArgument(4);
  unsigned base=0;
  unsigned nrows = myarg->getShape()[0];
  for(unsigned i=0; i<nrows; ++i) {
    unsigned ncols = myarg->getRowLength(i);
    for(unsigned j=0; j<ncols; ++j) {
      active_tasks.push_back(base+j);
    }
    base += myarg->getNumberOfColumns();
  }
  return active_tasks;
}

void QuaternionBondProductMatrix::performTask( const unsigned& taskno, MultiValue& myvals) const {
  unsigned index1 = std::floor( taskno / getPntrToArgument(4)->getNumberOfColumns() );
  unsigned index2 = taskno - getPntrToArgument(4)->getNumberOfColumns()*index1;

  std::vector<double> quat(4), bond(4), quatTemp(4);
  std::vector<Tensor4d> dqt(2); //dqt[0] -> derivs w.r.t quat [dwt/dw1 dwt/di1 dwt/dj1 dwt/dk1]
  //[dit/dw1 dit/di1 dit/dj1 dit/dk1] etc, and dqt[1] is w.r.t the vector-turned-quaternion called bond

  // Retrieve the quaternion
  for(unsigned i=0; i<4; ++i) {
    quat[i] = getPntrToArgument(i)->get(index1);
  }

  // Retrieve the components of the matrix
  double weight = getPntrToArgument(4)->get(taskno, false );
  for(unsigned i=1; i<4; ++i) {
    bond[i] = getPntrToArgument(4+i)->get(taskno, false );
  }

  // calculate normalization factor
  bond[0]=0.0;
  double normFac = 1/sqrt(bond[1]*bond[1] + bond[2]*bond[2] + bond[3]*bond[3]);
  if (bond[1] == 0.0 && bond[2]==0.0 && bond[3]==0) {
    normFac=1;  //just for the case where im comparing a quat to itself, itll be 0 at the end anyway
  }
  double normFac3 = normFac*normFac*normFac;
  //I hold off on normalizing because this can be done at the very end, and it makes the derivatives with respect to 'bond' more simple

  std::vector<double> quat_conj(4);
  quat_conj[0] = quat[0];
  quat_conj[1] = -1*quat[1];
  quat_conj[2] = -1*quat[2];
  quat_conj[3] = -1*quat[3];
  //make a conjugate of q1 my own sanity

  //q1_conj * r first, while keep track of derivs
  double pref=1;
  double conj=1;
  double pref2=1;
  //real part of q1*q2

  for(unsigned i=0; i<4; ++i) {
    if( i>0 ) {
      pref=-1;
      conj=-1;
      pref2=-1;
    }
    quatTemp[0]+=pref*quat_conj[i]*bond[i];
    dqt[0](0,i) = conj*pref*bond[i];
    dqt[1](0,i) = pref2*quat_conj[i];
  }
  //i component
  pref=1;
  conj=1;
  pref2=1;

  for (unsigned i=0; i<4; i++) {
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
    if (i>0) {
      conj=-1;
    }

    quatTemp[1]+=pref*quat_conj[i]*bond[(5-i)%4];
    dqt[0](1,i) =conj*pref*bond[(5-i)%4];
    dqt[1](1,i) = pref2*quat_conj[(5-i)%4];
  }

  //j component
  pref=1;
  pref2=1;
  conj=1;

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
    if (i>0) {
      conj=-1;
    }

    quatTemp[2]+=pref*quat_conj[i]*bond[(i+2)%4];
    dqt[0](2,i)=conj*pref*bond[(i+2)%4];
    dqt[1](2,i)=pref2*quat_conj[(i+2)%4];
  }

  //k component
  pref=1;
  pref2=1;
  conj=1;

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
    if(i>0) {
      conj=-1;
    }
    quatTemp[3]+=pref*quat_conj[i]*bond[(3-i)];
    dqt[0](3,i)=conj*pref*bond[3-i];
    dqt[1](3,i)= pref2*quat_conj[3-i];

  }


//now previous ^ product times quat again, not conjugated
  //real part of q1*q2
  double tempDot=0,wf=0,xf=0,yf=0,zf=0;
  pref=1;
  pref2=1;
  unsigned base=0;
  for(unsigned i=0; i<4; ++i) {
    if( i>0 ) {
      pref=-1;
      pref2=-1;
    }
    myvals.addValue( 0, normFac*pref*quatTemp[i]*quat[i] );
    wf+=normFac*pref*quatTemp[i]*quat[i];
    if( doNotCalculateDerivatives() ) {
      continue ;
    }
    tempDot=(dotProduct(Vector4d(quat[0],-quat[1],-quat[2],-quat[3]), dqt[0].getCol(i)) + pref2*quatTemp[i])*normFac;
    myvals.addDerivative( 0, base + index1, tempDot );
    myvals.updateIndex( 0, base + index1 );
    base += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  //had to split because bond's derivatives depend on the value of the overall quaternion component
  if( !doNotCalculateDerivatives() ) {
    for(unsigned i=0; i<4; ++i) {
      tempDot=dotProduct(Vector4d(quat[0],-quat[1],-quat[2],-quat[3]), dqt[1].getCol(i))*normFac;
      if( i>0 ) {
        myvals.addDerivative( 0, base + taskno, tempDot );
        myvals.updateIndex( 0, base + taskno );
      }
      base += getPntrToArgument(4+i)->getNumberOfStoredValues();
    }
  }

  //i component
  pref=1;
  pref2=1;
  base = 0;
  for (unsigned i=0; i<4; i++) {
    if(i==3) {
      pref=-1;
    } else {
      pref=1;
    }
    myvals.addValue( 1, normFac*pref*quatTemp[i]*quat[(5-i)%4]);
    xf+=normFac*pref*quatTemp[i]*quat[(5-i)%4];
    if(i==2) {
      pref2=-1;
    } else {
      pref2=1;
    }
    if( doNotCalculateDerivatives() ) {
      continue ;
    }
    tempDot=(dotProduct(Vector4d(quat[1],quat[0],quat[3],-quat[2]), dqt[0].getCol(i)) + pref2*quatTemp[(5-i)%4])*normFac;
    myvals.addDerivative( 1, base + index1, tempDot );
    myvals.updateIndex( 1, base + index1 );
    base += getPntrToArgument(i)->getNumberOfStoredValues();
  }

  if( !doNotCalculateDerivatives() ) {
    for(unsigned i=0; i<4; ++i) {
      tempDot=dotProduct(Vector4d(quat[1],quat[0],quat[3],-quat[2]), dqt[1].getCol(i))*normFac;
      if( i>0 ) {
        myvals.addDerivative( 1, base + taskno, tempDot+(-bond[i]*normFac*normFac*xf) );
        myvals.updateIndex( 1, base + taskno);
      }
      base += getPntrToArgument(4+i)->getNumberOfStoredValues();
    }
  }


  //j component
  pref=1;
  pref2=1;
  base = 0;
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

    myvals.addValue( 2, normFac*pref*quatTemp[i]*quat[(i+2)%4]);
    yf+=normFac*pref*quatTemp[i]*quat[(i+2)%4];
    if( doNotCalculateDerivatives() ) {
      continue ;
    }
    tempDot=(dotProduct(Vector4d(quat[2],-quat[3],quat[0],quat[1]), dqt[0].getCol(i)) + pref2*quatTemp[(i+2)%4])*normFac;
    myvals.addDerivative( 2, base + index1, tempDot );
    myvals.updateIndex( 2, base + index1 );
    base += getPntrToArgument(i)->getNumberOfStoredValues();
  }

  if( !doNotCalculateDerivatives() ) {
    for(unsigned i=0; i<4; ++i) {
      tempDot=dotProduct(Vector4d(quat[2],-quat[3],quat[0],quat[1]), dqt[1].getCol(i))*normFac;
      if( i>0 ) {
        myvals.addDerivative( 2, base + taskno, tempDot+(-bond[i]*normFac*normFac*yf) );
        myvals.updateIndex( 2, base + taskno );
      }
      base += getPntrToArgument(4+i)->getNumberOfStoredValues();
    }
  }

  //k component
  pref=1;
  pref2=1;
  base = 0;
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

    myvals.addValue( 3, normFac*pref*quatTemp[i]*quat[(3-i)]);
    zf+=normFac*pref*quatTemp[i]*quat[(3-i)];
    if( doNotCalculateDerivatives() ) {
      continue ;
    }
    tempDot=(dotProduct(Vector4d(quat[3],quat[2],-quat[1],quat[0]), dqt[0].getCol(i)) + pref2*quatTemp[(3-i)])*normFac;
    myvals.addDerivative( 3, base + index1, tempDot );
    myvals.updateIndex( 3, base + index1 );
    base += getPntrToArgument(i)->getNumberOfStoredValues();
  }

  if( doNotCalculateDerivatives() ) {
    return ;
  }

  for(unsigned i=0; i<4; ++i) {
    tempDot=dotProduct(Vector4d(quat[3],quat[2],-quat[1],quat[0]), dqt[1].getCol(i))*normFac;
    if( i>0 ) {
      myvals.addDerivative( 3, base + taskno, tempDot+(-bond[i]*normFac*normFac*zf) );
      myvals.updateIndex( 3, base + taskno);
    }
    base += getPntrToArgument(4+i)->getNumberOfStoredValues();
  }
}

void QuaternionBondProductMatrix::calculate() {
  // Copy bookeeping arrays from input matrices to output matrices
  for(unsigned i=0; i<4; ++i) {
    getPntrToComponent(i)->copyBookeepingArrayFromArgument( getPntrToArgument(4+i) );
  }
  runAllTasks();
}

}
}
