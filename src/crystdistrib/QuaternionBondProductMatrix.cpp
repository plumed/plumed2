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
#include "core/ActionWithMatrix.h"
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

class QuaternionBondProductMatrix : public ActionWithMatrix {
private:
  unsigned nderivatives;
  std::vector<bool> stored;
//  const Vector4d& rightMultiply(Tensor4d&, Vector4d&);
public:
  static void registerKeywords( Keywords& keys );
  explicit QuaternionBondProductMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  unsigned getNumberOfColumns() const override ;
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override;
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override ;
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
  ActionWithMatrix::registerKeywords(keys); keys.use("ARG");
  keys.addOutputComponent("w","default","the real component of quaternion");
  keys.addOutputComponent("i","default","the i component of the quaternion");
  keys.addOutputComponent("j","default","the j component of the quaternion");
  keys.addOutputComponent("k","default","the k component of the quaternion");
}

QuaternionBondProductMatrix::QuaternionBondProductMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao)
{
  if( getNumberOfArguments()!=8 ) error("should be eight arguments to this action, 4 quaternion components and 4 matrices");
  unsigned nquat = getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=0; i<4; ++i) {
    Value* myarg=getPntrToArgument(i); myarg->buildDataStore();
    if( myarg->getRank()!=1 ) error("first four arguments to this action should be vectors");
    if( (myarg->getPntrToAction())->getName()!="QUATERNION_VECTOR" ) error("first four arguments to this action should be quaternions");
    std::string mylab=getPntrToArgument(i)->getName(); std::size_t dot=mylab.find_first_of(".");
    if( i==0 && mylab.substr(dot+1)!="w" ) error("quaternion arguments are in wrong order");
    if( i==1 && mylab.substr(dot+1)!="i" ) error("quaternion arguments are in wrong order");
    if( i==2 && mylab.substr(dot+1)!="j" ) error("quaternion arguments are in wrong order");
    if( i==3 && mylab.substr(dot+1)!="k" ) error("quaternion arguments are in wrong order");
  }
  std::vector<unsigned> shape( getPntrToArgument(4)->getShape() );
  for(unsigned i=4; i<8; ++i) {
    Value* myarg=getPntrToArgument(i);
    if( myarg->getRank()!=2 ) error("second four arguments to this action should be matrices");
    if( myarg->getShape()[0]!=shape[0] || myarg->getShape()[1]!=shape[1] ) error("matrices should all have the same shape");
    if( myarg->getShape()[0]!=nquat ) error("number of rows in matrix should equal number of input quaternions");
    std::string mylab=getPntrToArgument(i)->getName(); std::size_t dot=mylab.find_first_of(".");
    if( i==5 && mylab.substr(dot+1)!="x" ) error("quaternion arguments are in wrong order");
    if( i==6 && mylab.substr(dot+1)!="y" ) error("quaternion arguments are in wrong order");
    if( i==7 && mylab.substr(dot+1)!="z" ) error("quaternion arguments are in wrong order");
  }
  addComponent( "w", shape ); componentIsNotPeriodic("w");
  addComponent( "i", shape ); componentIsNotPeriodic("i");
  addComponent( "j", shape ); componentIsNotPeriodic("j");
  addComponent( "k", shape ); componentIsNotPeriodic("k");
  done_in_chain=true; nderivatives = buildArgumentStore(0);

  std::string headstr=getFirstActionInChain()->getLabel(); stored.resize( getNumberOfArguments() );
  for(unsigned i=0; i<getNumberOfArguments(); ++i) stored[i] = getPntrToArgument(i)->ignoreStoredValue( headstr );
}

unsigned QuaternionBondProductMatrix::getNumberOfDerivatives() {
  return nderivatives;
}

unsigned QuaternionBondProductMatrix::getNumberOfColumns() const {
  const ActionWithMatrix* am=dynamic_cast<const ActionWithMatrix*>( getPntrToArgument(4)->getPntrToAction() );
  plumed_assert( am ); return am->getNumberOfColumns();
}

void QuaternionBondProductMatrix::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  unsigned start_n = getPntrToArgument(4)->getShape()[0], size_v = getPntrToArgument(4)->getShape()[1];
  if( indices.size()!=size_v+1 ) indices.resize( size_v+1 );
  for(unsigned i=0; i<size_v; ++i) indices[i+1] = start_n + i;
  myvals.setSplitIndex( size_v + 1 );
}

void QuaternionBondProductMatrix::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned ind2=index2;
  if( index2>=getPntrToArgument(0)->getShape()[0] ) ind2 = index2 - getPntrToArgument(0)->getShape()[0];

  std::vector<double> quat(4), bond(4), quatTemp(4);
  std::vector<Tensor4d> dqt(2); //dqt[0] -> derivs w.r.t quat [dwt/dw1 dwt/di1 dwt/dj1 dwt/dk1]
  //[dit/dw1 dit/di1 dit/dj1 dit/dk1] etc, and dqt[1] is w.r.t the vector-turned-quaternion called bond

  // Retrieve the quaternion
  for(unsigned i=0; i<4; ++i) quat[i] = getArgumentElement( i, index1, myvals );

  // Retrieve the components of the matrix
  double weight = getElementOfMatrixArgument( 4, index1, ind2, myvals );
  for(unsigned i=1; i<4; ++i) bond[i] = getElementOfMatrixArgument( 4+i, index1, ind2, myvals );

  // calculate normalization factor
  bond[0]=0.0;
  double normFac = 1/sqrt(bond[1]*bond[1] + bond[2]*bond[2] + bond[3]*bond[3]);
  if (bond[1] == 0.0 && bond[2]==0.0 && bond[3]==0) normFac=1; //just for the case where im comparing a quat to itself, itll be 0 at the end anyway
  double normFac3 = normFac*normFac*normFac;
  //I hold off on normalizing because this can be done at the very end, and it makes the derivatives with respect to 'bond' more simple



  std::vector<double> quat_conj(4);
  quat_conj[0] = quat[0]; quat_conj[1] = -1*quat[1]; quat_conj[2] = -1*quat[2]; quat_conj[3] = -1*quat[3];
  //make a conjugate of q1 my own sanity




//q1_conj * r first, while keep track of derivs
  double pref=1;
  double conj=1;
  double pref2=1;
  //real part of q1*q2

  for(unsigned i=0; i<4; ++i) {
    if( i>0 ) {pref=-1; conj=-1; pref2=-1;}
    quatTemp[0]+=pref*quat_conj[i]*bond[i];
    dqt[0](0,i) = conj*pref*bond[i];
    dqt[1](0,i) = pref2*quat_conj[i];
    //addDerivativeOnVectorArgument( false, 0, i, index1, conj*pref*bond[i], myvals );
    //addDerivativeOnVectorArgument( false, 0, 4+i, ind2, conj*pref*quat[i], myvals );
  }
  //i component
  pref=1;
  conj=1;
  pref2=1;

  for (unsigned i=0; i<4; i++) {
    if(i==3) pref=-1;
    else pref=1;
    if(i==2) pref2=-1;
    else pref2=1;
    if (i>0) conj=-1;

    quatTemp[1]+=pref*quat_conj[i]*bond[(5-i)%4];
    dqt[0](1,i) =conj*pref*bond[(5-i)%4];
    dqt[1](1,i) = pref2*quat_conj[(5-i)%4];
    //addDerivativeOnVectorArgument( false, 1, i, index1, conj*pref*bond[(5-i)%4], myvals );
    //addDerivativeOnVectorArgument( false, 1, 4+i, ind2, conj*pref*quat[i], myvals );
  }

  //j component
  pref=1;
  pref2=1;
  conj=1;

  for (unsigned i=0; i<4; i++) {
    if(i==1) pref=-1;
    else pref=1;
    if (i==3) pref2=-1;
    else pref2=1;
    if (i>0) conj=-1;

    quatTemp[2]+=pref*quat_conj[i]*bond[(i+2)%4];
    dqt[0](2,i)=conj*pref*bond[(i+2)%4];
    dqt[1](2,i)=pref2*quat_conj[(i+2)%4];
    //addDerivativeOnVectorArgument( false, 2, i, index1, conj*pref*bond[(i+2)%4], myvals );
    //addDerivativeOnVectorArgument( false, 2, 4+i, ind2, conj*pref*quat[i], myvals );
  }

  //k component
  pref=1;
  pref2=1;
  conj=1;

  for (unsigned i=0; i<4; i++) {
    if(i==2) pref=-1;
    else pref=1;
    if(i==1) pref2=-1;
    else pref2=1;
    if(i>0) conj=-1;
    quatTemp[3]+=pref*quat_conj[i]*bond[(3-i)];
    dqt[0](3,i)=conj*pref*bond[3-i];
    dqt[1](3,i)= pref2*quat_conj[3-i];
    //addDerivativeOnVectorArgument( false, 3, i, index1, conj*pref*bond[3-i], myvals );
    //addDerivativeOnVectorArgument( false, 3, 4+i, ind2, conj*pref*quat[i], myvals );

  }


//now previous ^ product times quat again, not conjugated
  //real part of q1*q2
  double tempDot=0,wf=0,xf=0,yf=0,zf=0;
  pref=1;
  pref2=1;
  for(unsigned i=0; i<4; ++i) {
    if( i>0 ) {pref=-1; pref2=-1;}
    myvals.addValue( getConstPntrToComponent(0)->getPositionInStream(), normFac*pref*quatTemp[i]*quat[i] );
    wf+=normFac*pref*quatTemp[i]*quat[i];
    if( doNotCalculateDerivatives() ) continue ;
    tempDot=(dotProduct(Vector4d(quat[0],-quat[1],-quat[2],-quat[3]), dqt[0].getCol(i)) + pref2*quatTemp[i])*normFac;
    addDerivativeOnVectorArgument( stored[i], 0, i,   index1, tempDot, myvals);
  }
  //had to split because bond's derivatives depend on the value of the overall quaternion component
  //addDerivativeOnMatrixArgument( false, 0, 4, index1, ind2, 0.0, myvals );
  for(unsigned i=0; i<4; ++i) {
    tempDot=dotProduct(Vector4d(quat[0],-quat[1],-quat[2],-quat[3]), dqt[1].getCol(i))*normFac;
    if (i!=0 )addDerivativeOnMatrixArgument( stored[4+i], 0, 4+i, index1, ind2, tempDot, myvals );
    else addDerivativeOnMatrixArgument( stored[4+i], 0, 4+i, index1, ind2, 0.0, myvals );
  }
// for (unsigned i=0; i<4; ++i) {
//myvals.addValue( getConstPntrToComponent(0)->getPositionInStream(), 0.0 );
//if( doNotCalculateDerivatives() ) continue ;
//addDerivativeOnVectorArgument( false, 0, i,   index1, 0.0, myvals);
//addDerivativeOnVectorArgument( false, 0, 4+i, ind2, 0.0 ,  myvals);
//  }
//the w component should always be zero, barring some catastrophe, but we calculate it out anyway

  //i component
  pref=1;
  pref2=1;
  for (unsigned i=0; i<4; i++) {
    if(i==3) pref=-1;
    else pref=1;
    myvals.addValue( getConstPntrToComponent(1)->getPositionInStream(), normFac*pref*quatTemp[i]*quat[(5-i)%4]);
    xf+=normFac*pref*quatTemp[i]*quat[(5-i)%4];
    if(i==2) pref2=-1;
    else pref2=1;
    if( doNotCalculateDerivatives() ) continue ;
    tempDot=(dotProduct(Vector4d(quat[1],quat[0],quat[3],-quat[2]), dqt[0].getCol(i)) + pref2*quatTemp[(5-i)%4])*normFac;
    addDerivativeOnVectorArgument( stored[i], 1, i,   index1, tempDot, myvals);
  }
  //addDerivativeOnMatrixArgument( false, 1, 4, index1, ind2, 0.0, myvals );

  for(unsigned i=0; i<4; ++i) {
    tempDot=dotProduct(Vector4d(quat[1],quat[0],quat[3],-quat[2]), dqt[1].getCol(i))*normFac;
    if (i!=0) addDerivativeOnMatrixArgument( stored[4+i], 1, 4+i, index1, ind2, tempDot+(-bond[i]*normFac*normFac*xf), myvals );
    else  addDerivativeOnMatrixArgument( stored[4+i], 1, 4+i, index1, ind2, 0.0, myvals );

  }


  //j component
  pref=1;
  pref2=1;
  for (unsigned i=0; i<4; i++) {
    if(i==1) pref=-1;
    else pref=1;
    if (i==3) pref2=-1;
    else pref2=1;

    myvals.addValue( getConstPntrToComponent(2)->getPositionInStream(), normFac*pref*quatTemp[i]*quat[(i+2)%4]);
    yf+=normFac*pref*quatTemp[i]*quat[(i+2)%4];
    if( doNotCalculateDerivatives() ) continue ;
    tempDot=(dotProduct(Vector4d(quat[2],-quat[3],quat[0],quat[1]), dqt[0].getCol(i)) + pref2*quatTemp[(i+2)%4])*normFac;
    addDerivativeOnVectorArgument( stored[i], 2, i,   index1, tempDot, myvals);
  }
  //    addDerivativeOnMatrixArgument( false, 2, 4, index1, ind2,0.0   , myvals );

  for(unsigned i=0; i<4; ++i) {
    tempDot=dotProduct(Vector4d(quat[2],-quat[3],quat[0],quat[1]), dqt[1].getCol(i))*normFac;
    if (i!=0) addDerivativeOnMatrixArgument( stored[4+i], 2, 4+i, index1, ind2, tempDot+(-bond[i]*normFac*normFac*yf), myvals );
    else  addDerivativeOnMatrixArgument( stored[4+i], 2, 4+i, index1, ind2, 0.0, myvals );


  }

  //k component
  pref=1;
  pref2=1;
  for (unsigned i=0; i<4; i++) {
    if(i==2) pref=-1;
    else pref=1;
    if(i==1) pref2=-1;
    else pref2=1;

    myvals.addValue( getConstPntrToComponent(3)->getPositionInStream(), normFac*pref*quatTemp[i]*quat[(3-i)]);
    zf+=normFac*pref*quatTemp[i]*quat[(3-i)];
    if( doNotCalculateDerivatives() ) continue ;
    tempDot=(dotProduct(Vector4d(quat[3],quat[2],-quat[1],quat[0]), dqt[0].getCol(i)) + pref2*quatTemp[(3-i)])*normFac;
    addDerivativeOnVectorArgument( stored[i], 3, i,   index1, tempDot, myvals);
  }
  //addDerivativeOnMatrixArgument( false, 3, 4, index1, ind2,  0.0 , myvals );

  for(unsigned i=0; i<4; ++i) {
    tempDot=dotProduct(Vector4d(quat[3],quat[2],-quat[1],quat[0]), dqt[1].getCol(i))*normFac;
    if (i!=0) addDerivativeOnMatrixArgument( stored[4+i], 3, 4+i, index1, ind2, tempDot+(-bond[i]*normFac*normFac*zf), myvals );
    else addDerivativeOnMatrixArgument( stored[4+i], 3, 4+i, index1, ind2, 0.0, myvals );


  }
  if( doNotCalculateDerivatives() ) return ;

  for(unsigned outcomp=0; outcomp<4; ++outcomp) {
    unsigned ostrn = getConstPntrToComponent(outcomp)->getPositionInStream();
    for(unsigned i=4; i<8; ++i) {
      bool found=false;
      for(unsigned j=4; j<i; ++j) {
        if( arg_deriv_starts[i]==arg_deriv_starts[j] ) { found=true; break; }
      }
      if( found || !stored[i] ) continue;

      unsigned istrn = getPntrToArgument(i)->getPositionInStream();
      for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
        unsigned kind=myvals.getActiveIndex(istrn,k); myvals.updateIndex( ostrn, kind );
      }
    }
  }
}

void QuaternionBondProductMatrix::runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() || !matrixChainContinues() ) return ;

  for(unsigned j=0; j<getNumberOfComponents(); ++j) {
    unsigned nmat = getConstPntrToComponent(j)->getPositionInMatrixStash(), nmat_ind = myvals.getNumberOfMatrixRowDerivatives( nmat );
    std::vector<unsigned>& matrix_indices( myvals.getMatrixRowDerivativeIndices( nmat ) ); unsigned ntwo_atoms = myvals.getSplitIndex();
    // Quaternion
    for(unsigned k=0; k<4; ++k) { matrix_indices[nmat_ind] = arg_deriv_starts[k] + ival; nmat_ind++; }
    // Loop over row of matrix
    for(unsigned n=4; n<8; ++n) {
      bool found=false;
      for(unsigned k=4; k<n; ++k) {
        if( arg_deriv_starts[k]==arg_deriv_starts[n] ) { found=true; break; }
      }
      if( found ) continue;
      unsigned istrn = getPntrToArgument(n)->getPositionInMatrixStash();
      std::vector<unsigned>& imat_indices( myvals.getMatrixRowDerivativeIndices( istrn ) );
      for(unsigned k=0; k<myvals.getNumberOfMatrixRowDerivatives( istrn ); ++k) matrix_indices[nmat_ind + k] = arg_deriv_starts[n] + imat_indices[k];
      nmat_ind += myvals.getNumberOfMatrixRowDerivatives( getPntrToArgument(4)->getPositionInMatrixStash() );
    }
    myvals.setNumberOfMatrixRowDerivatives( nmat, nmat_ind );
  }
}

}
}
