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
#include "core/ParallelTaskManager.h"
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

class QuatBondProdMatInput {
public:
  int fake{1};
};

class QuaternionBondProductMatrix : public ActionWithVector {
public:
  using input_type = QuatBondProdMatInput;
  using PTM = ParallelTaskManager<QuaternionBondProductMatrix>;
private:
  PTM taskmanager;
  std::vector<unsigned> active_tasks;
//  const Vector4d& rightMultiply(Tensor4d&, Vector4d&);
public:
  static void registerKeywords( Keywords& keys );
  explicit QuaternionBondProductMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  void prepare() override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  std::vector<unsigned>& getListOfActiveTasks( ActionWithVector* action ) override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const override {
    plumed_merror("this doesn't do anything");
  }
  static void performTask( std::size_t task_index,
                           const QuatBondProdMatInput& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  static void gatherForces( std::size_t task_index,
                            const QuatBondProdMatInput& actiondata,
                            const ParallelActionsInput& input,
                            const ForceInput& fdata,
                            ForceOutput forces );
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
  PTM::registerKeywords( keys );
  keys.addOutputComponent("w","default","matrix","the real component of quaternion");
  keys.addOutputComponent("i","default","matrix","the i component of the quaternion");
  keys.addOutputComponent("j","default","matrix","the j component of the quaternion");
  keys.addOutputComponent("k","default","matrix","the k component of the quaternion");
}

QuaternionBondProductMatrix::QuaternionBondProductMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this) {
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
  std::vector<std::size_t> shape( getPntrToArgument(4)->getShape() );
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
  // We now swap round the order of the arguments to make gathering of forces easier for parallel task manager
  std::vector<Value*> args( getArguments() ), newargs;
  for(unsigned i=0; i<4; ++i) {
    newargs.push_back(args[4+i]);
  }
  for(unsigned i=0; i<4; ++i) {
    newargs.push_back(args[i]);
  }
  requestArguments( newargs );

  addComponent( "w", shape );
  componentIsNotPeriodic("w");
  addComponent( "i", shape );
  componentIsNotPeriodic("i");
  addComponent( "j", shape );
  componentIsNotPeriodic("j");
  addComponent( "k", shape );
  componentIsNotPeriodic("k");
  taskmanager.setupParallelTaskManager( 1, 8, 4*nquat );
  taskmanager.setActionInput( QuatBondProdMatInput() );
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

  Value* myarg = getPntrToArgument(0);
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

void QuaternionBondProductMatrix::performTask( std::size_t task_index,
    const QuatBondProdMatInput& actiondata,
    ParallelActionsInput& input,
    ParallelActionsOutput& output ) {

  View2D<double, 4, 8> derivatives( output.derivatives.data() );
  std::size_t index1 = std::floor( task_index / input.args[0].ncols );
  std::vector<double> quat(4), bond(4), quatTemp(4);
  std::vector<Tensor4d> dqt(2); //dqt[0] -> derivs w.r.t quat [dwt/dw1 dwt/di1 dwt/dj1 dwt/dk1]
  //[dit/dw1 dit/di1 dit/dj1 dit/dk1] etc, and dqt[1] is w.r.t the vector-turned-quaternion called bond

  // Retrieve the quaternion
  for(unsigned i=0; i<4; ++i) {
    quat[i] = input.inputdata[ input.args[4+i].start + index1 ];
  }

  // Retrieve the components of the matrix
  double weight = input.inputdata[ input.args[0].start + task_index ];
  for(unsigned i=1; i<4; ++i) {
    bond[i] = input.inputdata[ input.args[i].start + task_index ];
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
  for(unsigned i=0; i<4; ++i) {
    if( i>0 ) {
      pref=-1;
      pref2=-1;
    }
    output.values[0] += normFac*pref*quatTemp[i]*quat[i];
    wf+=normFac*pref*quatTemp[i]*quat[i];
    if( input.noderiv ) {
      continue ;
    }
    tempDot=(dotProduct(Vector4d(quat[0],-quat[1],-quat[2],-quat[3]), dqt[0].getCol(i)) + pref2*quatTemp[i])*normFac;
    derivatives[0][4+i] = tempDot;
  }
  //had to split because bond's derivatives depend on the value of the overall quaternion component
  if( !input.noderiv ) {
    for(unsigned i=0; i<4; ++i) {
      tempDot=dotProduct(Vector4d(quat[0],-quat[1],-quat[2],-quat[3]), dqt[1].getCol(i))*normFac;
      if( i>0 ) {
        derivatives[0][i] = tempDot;
      }
    }
  }

  //i component
  pref=1;
  pref2=1;
  for (unsigned i=0; i<4; i++) {
    if(i==3) {
      pref=-1;
    } else {
      pref=1;
    }
    output.values[1] += normFac*pref*quatTemp[i]*quat[(5-i)%4];
    xf+=normFac*pref*quatTemp[i]*quat[(5-i)%4];
    if(i==2) {
      pref2=-1;
    } else {
      pref2=1;
    }
    if( input.noderiv ) {
      continue ;
    }
    tempDot=(dotProduct(Vector4d(quat[1],quat[0],quat[3],-quat[2]), dqt[0].getCol(i)) + pref2*quatTemp[(5-i)%4])*normFac;
    derivatives[1][4+i] = tempDot;
  }

  if( !input.noderiv ) {
    for(unsigned i=0; i<4; ++i) {
      tempDot=dotProduct(Vector4d(quat[1],quat[0],quat[3],-quat[2]), dqt[1].getCol(i))*normFac;
      if( i>0 ) {
        derivatives[1][i] = tempDot+(-bond[i]*normFac*normFac*xf);
      }
    }
  }

  //j component
  pref=1;
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

    output.values[2] += normFac*pref*quatTemp[i]*quat[(i+2)%4];
    yf+=normFac*pref*quatTemp[i]*quat[(i+2)%4];
    if( input.noderiv ) {
      continue ;
    }
    tempDot=(dotProduct(Vector4d(quat[2],-quat[3],quat[0],quat[1]), dqt[0].getCol(i)) + pref2*quatTemp[(i+2)%4])*normFac;
    derivatives[2][4+i] = tempDot;
  }

  if( !input.noderiv ) {
    for(unsigned i=0; i<4; ++i) {
      tempDot=dotProduct(Vector4d(quat[2],-quat[3],quat[0],quat[1]), dqt[1].getCol(i))*normFac;
      if( i>0 ) {
        derivatives[2][i] = tempDot+(-bond[i]*normFac*normFac*yf);
      }
    }
  }

//k component
  pref=1;
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

    output.values[3] += normFac*pref*quatTemp[i]*quat[(3-i)];
    zf+=normFac*pref*quatTemp[i]*quat[(3-i)];
    if( input.noderiv ) {
      continue ;
    }
    tempDot=(dotProduct(Vector4d(quat[3],quat[2],-quat[1],quat[0]), dqt[0].getCol(i)) + pref2*quatTemp[(3-i)])*normFac;
    derivatives[3][4+i] = tempDot;
  }

  if( input.noderiv ) {
    return ;
  }

  for(unsigned i=0; i<4; ++i) {
    tempDot=dotProduct(Vector4d(quat[3],quat[2],-quat[1],quat[0]), dqt[1].getCol(i))*normFac;
    if( i>0 ) {
      derivatives[3][i] = tempDot+(-bond[i]*normFac*normFac*zf);
    }
  }
}

void QuaternionBondProductMatrix::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

void QuaternionBondProductMatrix:: gatherForces( std::size_t task_index,
    const QuatBondProdMatInput& actiondata,
    const ParallelActionsInput& input,
    const ForceInput& fdata,
    ForceOutput forces ) {
  std::size_t nquat = input.args[4].shape[0];
  std::size_t index1 = std::floor( task_index / input.args[0].ncols );
  for(unsigned i=0; i<4; ++i) {
    double ff = fdata.force[i];
    for(unsigned j=0; j<4; ++j) {
      forces.thread_unsafe[ input.args[j].start + task_index] += ff*fdata.deriv[i][j];
    }
    for(unsigned j=0; j<4; ++j) {
      forces.thread_safe[ j*nquat + index1 ] += ff*fdata.deriv[i][4+j];
    }
  }
}

void QuaternionBondProductMatrix::calculate() {
  // Copy bookeeping arrays from input matrices to output matrices
  for(unsigned i=0; i<4; ++i) {
    getPntrToComponent(i)->copyBookeepingArrayFromArgument( getPntrToArgument(i) );
  }
  taskmanager.runAllTasks();
}

}
}
