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
// #ifdef __PLUMED_HAS_OPENACC
// #define __PLUMED_USE_OPENACC 1
// #endif //__PLUMED_HAS_OPENACC
#include "core/ActionWithVector.h"
#include "core/ParallelTaskManager.h"
#include "core/ActionRegister.h"
#include "tools/Torsion.h"


#include <iostream>

namespace PLMD {
namespace crystdistrib {

//+PLUMEDOC MCOLVAR QUATERNION_BOND_PRODUCT_MATRIX
/*
Calculate the product between a matrix of quaternions and the bonds connecting molecules

Calculate the product of a quaternion, times a vector, times the conjugate of the quaternion. Geometrically, this is the given vector projected onto the coordinate system defined by the quaternion. For context, the QUATERNION action defines an orthogonal coordinate frame given 3 atoms (i.e. one can define a coordinate system based on the structure of a given molecule). This action can then project a vector from the given molecular frame, toward another molecule, essentially pointing toward molecule 2, from the perspective of molecule 1. See QUATERNION for information about the molecular coordinate frame. Given a quaternion centered on molecule 1 $\mathbf{q1}$, and the vector connecting molecule 1, and some other molecule 2, $\mathbf{r_{21}}$, the following calculation is performed:

$$
\mathbf{r} = \overline{\mathbf{q_1}} \mathbf{r_{21}} \mathbf{q_1}
$$

where the overline denotes the quaternion conjugate. Internally, the vector $\mathbf{r_{21}}$ is treated as a quaternion with zero real part. Such a multiplication will always yield another quaternion with zero real part, and the results can be interpreted as an ordinary vector in $\mathbb{R}^3$ Nevertheless, this action has four components, the first of which, w, will always be entirely zeros. Finally, the resulting vector is normalized within the action, and the real length can be returned by multiplying each component by the norm of the vector given to the action. The quaternion should be a vector, and the distances a matrix.

In this example, 3 quaternion frames are calculated, and multiplied element-wise onto a distance matrix, yielding 9 vectors overall.

```plumed
#calculate the quaternion frames for 3 molecules
quat: QUATERNION ATOMS1=1,2,3 ATOMS2=4,5,6 ATOMS3=7,8,9
#also find the distance between the 3 origins of the molecule frames
c1: DISTANCE_MATRIX GROUP=1,4,7 CUTOFF=100.0 COMPONENTS
qp: QUATERNION_BOND_PRODUCT_MATRIX ARG=quat.*,c1.*
#this is now a matrix showing how each molecule is oriented in 3D space
#relative to eachother's origins
#now use in a simple collective variable
ops: CUSTOM ARG=qp.* VAR=w,i,j,k FUNC=w+i+j+k PERIODIC=NO
#w could have been ignored because it is always zero.

```

*/
//+ENDPLUMEDOC

struct QuatBondProdMatInput {
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(this[0:1])
  }

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
  void getNumberOfTasks( unsigned& ntasks ) override ;
  std::vector<unsigned>& getListOfActiveTasks( ActionWithVector* action ) override ;
  static void performTask( std::size_t task_index,
                           const QuatBondProdMatInput& actiondata,
                           const ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  static int getNumberOfValuesPerTask( std::size_t task_index, const QuatBondProdMatInput& actiondata );
  static void getForceIndices( std::size_t task_index, std::size_t colno, std::size_t ntotal_force, const QuatBondProdMatInput& actiondata, const ParallelActionsInput& input, ForceIndexHolder force_indices );
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
  std::size_t nonthreadsafe = 0;
  for(unsigned i=4; i<8; ++i) {
    nonthreadsafe += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  taskmanager.setupParallelTaskManager( 8, nonthreadsafe );
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

void QuaternionBondProductMatrix::getNumberOfTasks( unsigned& ntasks ) {
  ntasks=getPntrToComponent(0)->getNumberOfStoredValues();
}

std::vector<unsigned>& QuaternionBondProductMatrix::getListOfActiveTasks( ActionWithVector* action ) {
  if( active_tasks.size()>0 ) {
    return active_tasks;
  }

  unsigned atsize = 0;
  Value* myarg = getPntrToArgument(0);
  unsigned nrows = myarg->getShape()[0];
  for(unsigned i=0; i<nrows; ++i) {
    atsize += myarg->getRowLength(i);
  }
  active_tasks.resize( atsize );

  unsigned base=0, k=0;
  for(unsigned i=0; i<nrows; ++i) {
    unsigned ncols = myarg->getRowLength(i);
    for(unsigned j=0; j<ncols; ++j) {
      active_tasks[k] = base+j;
      ++k;
    }
    base += myarg->getNumberOfColumns();
  }
  return active_tasks;
}

void QuaternionBondProductMatrix::performTask( std::size_t task_index,
    const QuatBondProdMatInput& /*actiondata*/,
    const ParallelActionsInput& input,
    ParallelActionsOutput& output ) {

  View2D<double, 4, 8> derivatives( output.derivatives.data() );
  std::size_t index1 = std::floor( task_index / input.ncols[0] );

  std::array<Tensor4d,2> dqt; //dqt[0] -> derivs w.r.t quat [dwt/dw1 dwt/di1 dwt/dj1 dwt/dk1]
  //[dit/dw1 dit/di1 dit/dj1 dit/dk1] etc, and dqt[1] is w.r.t the vector-turned-quaternion called bond

  // Retrieve the quaternion
  const std::array<double,4> quat{
    input.inputdata[ input.argstarts[4+0] + index1 ],
    input.inputdata[ input.argstarts[4+1] + index1 ],
    input.inputdata[ input.argstarts[4+2] + index1 ],
    input.inputdata[ input.argstarts[4+3] + index1 ]
  };
  // Retrieve the components of the matrix
  const std::array<double,4> bond{
    0.0,
    input.inputdata[ input.argstarts[1] + task_index ],
    input.inputdata[ input.argstarts[2] + task_index ],
    input.inputdata[ input.argstarts[3] + task_index ]
  };

  //make a conjugate of q1 my own sanity
  Vector4d quat_conj{ quat[0],-quat[1],-quat[2],-quat[3]};
  //declaring some signs to simplify the loop
  constexpr std::array<double,4> pref_i{1,1,1,-1};
  constexpr std::array<double,4> pref2_i{1,1,-1,1};

  constexpr std::array<double,4> pref_j{1,-1,1,1};
  constexpr std::array<double,4> pref2_j{1,1,1,-1};

  constexpr std::array<double,4> pref_k{1,1,-1,1};
  constexpr std::array<double,4> pref2_k{1,-1,1,1};

  constexpr std::array<double,4> conj{1,-1,-1,-1};
  std::array<double,4> quatTemp{0.0,0.0,0.0,0.0};
  //q1_conj * r first, while keep track of derivs
  for(unsigned i=0; i<4; ++i) {
    //real part of q1*q2
    quatTemp[0]+= quat[i]*bond[i];//conj*pref 1*1 for 0, -1 * -1 for 1++
    dqt[0](0,i) = bond[i];//conj*pref 1*1 for 0, -1 * -1 for 1++
    dqt[1](0,i) = quat[i];//conj*pref 1*1 for 0, -1 * -1 for 1++
    //i component
    quatTemp[1]+= pref_i[i]*quat_conj[i]*bond[(5-i)%4];
    dqt[0](1,i) = pref_i[i]*conj[i]*bond[(5-i)%4];
    dqt[1](1,i) = pref2_i[i]*quat_conj[(5-i)%4];
    //j component
    quatTemp[2]+= pref_j[i]*quat_conj[i]*bond[(i+2)%4];
    dqt[0](2,i) = pref_j[i]*conj[i]*bond[(i+2)%4];
    dqt[1](2,i) = pref2_j[i]*quat_conj[(i+2)%4];
    //k component
    quatTemp[3]+= pref_k[i]*quat_conj[i]*bond[3-i];
    dqt[0](3,i) = pref_k[i]*conj[i]*bond[3-i];
    dqt[1](3,i) = pref2_k[i]*quat_conj[3-i];
  }
//now previous ^ product times quat again, not conjugated

// calculate normalization factor
  const double normFac = (bond[1] == 0.0 && bond[2]==0.0 && bond[3]==0) ?
                         //just for the case where im comparing a quat to itself, itll be 0 at the end anyway
                         1: 1/sqrt(bond[1]*bond[1] + bond[2]*bond[2] + bond[3]*bond[3]);
//I hold off on normalizing because this can be done at the very end, and it
// makes the derivatives with respect to 'bond' more simple

  double wf=0,xf=0,yf=0,zf=0;

  for(unsigned i=0; i<4; ++i) {
    //real part of q1*q2
    output.values[0] += normFac*conj[i]*quatTemp[i]*quat[i];
    wf+=normFac*conj[i]*quatTemp[i]*quat[i];
    //i component
    output.values[1] += normFac*pref_i[i]*quatTemp[i]*quat[(5-i)%4];
    xf+=normFac*pref_i[i]*quatTemp[i]*quat[(5-i)%4];
    //j component
    output.values[2] += normFac*pref_j[i]*quatTemp[i]*quat[(i+2)%4];
    yf+=normFac*pref_j[i]*quatTemp[i]*quat[(i+2)%4];
    //k component
    output.values[3] += normFac*pref_k[i]*quatTemp[i]*quat[(3-i)];
    zf+=normFac*pref_k[i]*quatTemp[i]*quat[(3-i)];
  }
  //had to split because bond's derivatives depend on the value of the overall quaternion component
  if( input.noderiv ) {
    return ;
  }
  const auto normFacSQ=normFac*normFac;
  xf*=normFacSQ;
  yf*=normFacSQ;
  zf*=normFacSQ;
  {
    //I unrolled the loop to opt out the ifs
    //I copy pasted the first component and constexpr'd the index
    constexpr int i = 0;
    //real part of q1*q2
    derivatives[0][4+i] = (dotProduct(quat_conj, dqt[0].getCol(i)) + conj[i]*quatTemp[i])*normFac;
    //i component
    derivatives[1][4+i] = (dotProduct(Vector4d(quat[1],quat[0],quat[3],-quat[2]),
                                      dqt[0].getCol(i)) + pref2_i[i]*quatTemp[(5-i)%4])*normFac;
    //j component
    derivatives[2][4+i] = (dotProduct(Vector4d(quat[2],-quat[3],quat[0],quat[1]),
                                      dqt[0].getCol(i)) + pref2_j[i]*quatTemp[(i+2)%4])*normFac;
    //k component
    derivatives[3][4+i] = (dotProduct(Vector4d(quat[3],quat[2],-quat[1],quat[0]),
                                      dqt[0].getCol(i)) + pref2_k[i]*quatTemp[(3-i)])*normFac;
  }
  for(unsigned i=1; i<4; ++i) {
    //real part of q1*q2
    derivatives[0][4+i] = (dotProduct(quat_conj, dqt[0].getCol(i)) + conj[i]*quatTemp[i])*normFac;
    derivatives[0][i]   = dotProduct(quat_conj, dqt[1].getCol(i))*normFac;
    //i component
    derivatives[1][4+i] =(dotProduct(Vector4d(quat[1],quat[0],quat[3],-quat[2]),
                                     dqt[0].getCol(i)) + pref2_i[i]*quatTemp[(5-i)%4])*normFac;
    derivatives[1][i]   = dotProduct(Vector4d(quat[1],quat[0],quat[3],-quat[2]),
                                     dqt[1].getCol(i))*normFac
                          +(-bond[i]*xf);
    //j component
    derivatives[2][4+i] = (dotProduct(Vector4d(quat[2],-quat[3],quat[0],quat[1]),
                                      dqt[0].getCol(i)) + pref2_j[i]*quatTemp[(i+2)%4])*normFac;
    derivatives[2][i] = dotProduct(Vector4d(quat[2],-quat[3],quat[0],quat[1]),
                                   dqt[1].getCol(i))*normFac+(-bond[i]*yf);
    //k component
    derivatives[3][4+i] = (dotProduct(Vector4d(quat[3],quat[2],-quat[1],quat[0]),
                                      dqt[0].getCol(i)) + pref2_k[i]*quatTemp[(3-i)])*normFac;
    derivatives[3][i] = dotProduct(Vector4d(quat[3],quat[2],-quat[1],quat[0]),
                                   dqt[1].getCol(i))*normFac+(-bond[i]*zf);
  }
}

void QuaternionBondProductMatrix::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

int QuaternionBondProductMatrix::getNumberOfValuesPerTask( std::size_t task_index, const QuatBondProdMatInput& actiondata ) {
  return 1;
}

void QuaternionBondProductMatrix::getForceIndices( std::size_t task_index,
    std::size_t colno,
    std::size_t ntotal_force,
    const QuatBondProdMatInput& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {
  std::size_t nquat = input.shapedata[input.shapestarts[4]];
  std::size_t index1 = std::floor( task_index / input.ncols[0] );
  for(unsigned i=0; i<4; ++i) {
    for(unsigned j=0; j<4; ++j) {
      force_indices.indices[i][j] = input.argstarts[j] + task_index;
    }
    force_indices.threadsafe_derivatives_end[i] = 4;
    for(unsigned j=0; j<4; ++j) {
      force_indices.indices[i][j+4] = j*nquat + index1;
    }
    force_indices.tot_indices[i] = 8;
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
