/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "adjmat/MatrixProductBase.h"
#include "core/ActionRegister.h"
#include "tools/LeptonCall.h"
#include "tools/Angle.h"

namespace PLMD {
namespace symfunc {

class ThreeBodyGFunctions : public adjmat::MatrixProductBase {
private:
  unsigned nderivatives;
  std::vector<LeptonCall> functions;
public:
  static void registerKeywords( Keywords& keys );
  explicit ThreeBodyGFunctions(const ActionOptions&);
  unsigned getNumberOfDerivatives() const override;
  void performTask( const unsigned& task_index, MultiValue& myvals ) const override;
  bool performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override {};
  double computeVectorProduct( const unsigned& index1, const unsigned& index2,
                                       const std::vector<double>& vec1, const std::vector<double>& vec2,
                                       std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const override { plumed_error(); }
};

PLUMED_REGISTER_ACTION(ThreeBodyGFunctions,"GSYMFUNC_THREEBODY")

void ThreeBodyGFunctions::registerKeywords( Keywords& keys ) {
  adjmat::MatrixProductBase::registerKeywords( keys );
  keys.add("compulsory","WEIGHT","the matrix that contains the weights that should be used for each connection"); 
  keys.add("numbered","FUNCTION","the parameters of the function you would like to compute");
  keys.add("compulsory","SWITCH","the switching function that acts on the distance between atom j and atom k in the G4 symmetry function");
  ActionWithValue::useCustomisableComponents( keys );
}

ThreeBodyGFunctions::ThreeBodyGFunctions(const ActionOptions&ao):
  Action(ao),
  MatrixProductBase(ao)
{
  if( getNumberOfArguments()!=3 ) error("found wrong number of arguments in input");
  std::vector<Value*> wval; parseArgumentList("WEIGHT",wval); 
  if( wval.size()!=1 ) error("keyword WEIGHT should be provided with the label of a single action");
  nderivatives=(wval[0]->getPntrToAction())->getNumberOfDerivatives();

  for(unsigned i=0;i<3;++i) {
      if( getPntrToArgument(i)->getRank()!=2 ) error("input argument should be a matrix");
      if( wval[0]->getShape()[0]!=getPntrToArgument(i)->getShape()[0] || wval[0]->getShape()[1]!=getPntrToArgument(i)->getShape()[1] ) error("mismatched shapes of matrices in input");
      bool found=false; if( wval[0]->getPntrToAction()->getLabel()==getPntrToArgument(i)->getPntrToAction()->getLabel() ) found=true; 
      if (!found ) {
          for(unsigned j=0;j<i;++j) { 
              if( getPntrToArgument(i)->getPntrToAction()->getLabel()==getPntrToArgument(i)->getPntrToAction()->getLabel() ) { found=true; break; }
          }
          if( !found ) nderivatives += (getPntrToArgument(i)->getPntrToAction())->getNumberOfDerivatives(); 
      }
  }
  std::vector<std::string> alabels(1); alabels[0]=(wval[0]->getPntrToAction())->getLabel();
  (wval[0]->getPntrToAction())->addActionToChain( alabels, this );
  log.printf("  using bond weights from matrix labelled %s \n",wval[0]->getName().c_str() );  
  // Rerequest the arguments
  std::vector<Value*> myargs( getArguments() ); myargs.push_back( wval[0] ); requestArguments( myargs, true );
  for(unsigned i=0; i<myargs.size(); ++i) myargs[i]->buildDataStore( getLabel() + "_mat" );
  std::vector<unsigned> shape(1); shape[0] = getPntrToArgument(0)->getShape()[0];

  // And now read the functions to compute
  for(int i=1;; i++) {
    std::string myfunc, mystr, lab, num; Tools::convert(i,num);
    if( !parseNumbered("FUNCTION",i,mystr ) ) break;
    std::vector<std::string> data=Tools::getWords(mystr);
    if( !Tools::parse(data,"LABEL",lab ) ) error("found no LABEL in FUNCTION" + num + " specification");
    addComponent( lab, shape ); componentIsNotPeriodic( lab );
    if( !Tools::parse(data,"FUNC",myfunc) ) error("found no FUNC in FUNCTION" + num + " specification");
    log.printf("  component labelled %s is computed using %s \n",lab.c_str(), myfunc.c_str() ); 
    functions.push_back( LeptonCall() ); std::vector<std::string> argnames(1); argnames[0]="ajik";
    if( myfunc.find("rij")!=std::string::npos ) argnames.push_back("rij"); 
    if( myfunc.find("rik")!=std::string::npos ) { 
        if( argnames.size()<2 ) error("if you have a function of rik it must also be a function of rij -- email gareth.tribello@gmail.com if this is a problem");
        argnames.push_back("rik"); 
    }
    if( myfunc.find("rjk")!=std::string::npos ) {
        if( argnames.size()<2 ) error("if you have a function of rjk it must also be a function of rij and rik -- email gareth.tribello@gmail.com if this is a problem");
        argnames.push_back("rjk"); 
    }
    functions[i-1].set( myfunc, argnames, this, true ); 
  }
  checkRead();
}

unsigned ThreeBodyGFunctions::getNumberOfDerivatives() const {
  return nderivatives;
}

void ThreeBodyGFunctions::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  if( !myvals.inVectorCall() ) return ;

  std::vector<Vector> der_i(4), der_j(4); 
  std::vector<double> values(4); Angle angle; Vector disti, distj; 
  plumed_assert( getNumberOfArguments()==4 );
  unsigned matind = getPntrToArgument(3)->getPositionInMatrixStash();
  unsigned matind_x = getPntrToArgument(0)->getPositionInMatrixStash();
  unsigned matind_y = getPntrToArgument(1)->getPositionInMatrixStash();
  unsigned matind_z = getPntrToArgument(2)->getPositionInMatrixStash();

  std::vector<std::vector<double> > deriv( functions.size() ); 
  if( !doNotCalculateDerivatives() ) {
      for(unsigned i=0; i<functions.size(); ++i) deriv[i].resize( 4*myvals.getNumberOfStashedMatrixElements(matind), 0 );
  }

  for(unsigned i=1; i<myvals.getNumberOfStashedMatrixElements(matind); ++i) {
    unsigned iind = myvals.getStashedMatrixIndex(matind,i);
    double weighti = myvals.getStashedMatrixElement( matind, iind );
    if( weighti<epsilon ) continue ;
    disti[0] = myvals.getStashedMatrixElement( matind_x, iind );
    disti[1] = myvals.getStashedMatrixElement( matind_y, iind );
    disti[2] = myvals.getStashedMatrixElement( matind_z, iind );
    values[1] = disti.modulo2(); der_i[1]=2*disti; der_i[2].zero();
    for(unsigned j=0; j<i; ++j) {
      unsigned jind = myvals.getStashedMatrixIndex(matind,j);
      double weightj = myvals.getStashedMatrixElement( matind, jind );
      if( weightj<epsilon ) continue ;
      distj[0] = myvals.getStashedMatrixElement( matind_x, jind );
      distj[1] = myvals.getStashedMatrixElement( matind_y, jind );
      distj[2] = myvals.getStashedMatrixElement( matind_z, jind );
      values[2] = distj.modulo2(); der_j[1].zero(); der_j[2]=2*distj;
      der_i[3] = disti - distj; values[3] = der_i[3].modulo2(); 
      der_i[3] = 2*der_i[3]; der_j[3] = -der_i[3];
      // Compute angle between bonds
      values[0] = angle.compute( disti, distj, der_i[0], der_j[0] );
      // Compute product of weights
      double weightij = weighti*weightj;
      // Now compute all symmetry functions
      for(unsigned n=0; n<functions.size(); ++n) {
        unsigned ostrn = getPntrToOutput(n)->getPositionInStream();
        double nonweight = functions[n].evaluate( values ); myvals.addValue( ostrn, nonweight*weightij ); 
        if( doNotCalculateDerivatives() ) continue;

        for(unsigned m=0;m<functions[n].getNumberOfArguments();++m) {
            unsigned base=0; double der = weightij*functions[n].evaluateDeriv( m, values );
            for(unsigned k=0; k<3; ++k) {
                deriv[n][base+i] += der*der_i[m][k]; deriv[n][base+j] += der*der_j[m][k];
                base += myvals.getNumberOfStashedMatrixElements(matind);
            }
        }
        unsigned base=3*myvals.getNumberOfStashedMatrixElements(matind);
        deriv[n][base+i] += nonweight*weightj; deriv[n][base+j] += nonweight*weighti;
      }
    }
  }
  if( doNotCalculateDerivatives() ) return;
  // First turn off matrix element storing during rerun of calculation
  myvals.setMatrixStashForRerun();
  // Now calculate the derivatives by recomputing the matrix elements that have derivatives 
  unsigned aindex_start = myvals.getNumberOfIndicesInFirstBlock(); 
  unsigned nstash = myvals.getNumberOfStashedMatrixElements(matind);
  ActionWithValue* av = (getPntrToArgument(3)->getPntrToAction())->getActionThatCalculates();
  for(unsigned i=0; i<nstash; ++i) {
      // Check there are derivatives
      double totder = 0;
      for(unsigned n=0; n<functions.size(); ++n) {
          for(unsigned k=0; k<4; ++k) totder += deriv[n][ k*nstash + i ];
      }
      if( fabs(totder)<epsilon ) continue;

      unsigned iind = myvals.getStashedMatrixIndex(matind,i);
      // Rerun the task
      av->runTask( av->getLabel(), task_index, aindex_start + iind, myvals );
      // Add the derivatives
      for(unsigned n=0; n<functions.size(); ++n) {
          unsigned ostrn = getPntrToOutput(n)->getPositionInStream();
          for(unsigned k=0; k<4; ++k) {
              double der = deriv[n][ k*nstash + i];
              unsigned istrn = getPntrToArgument(k)->getPositionInStream(); 
              for(unsigned j=0; j<myvals.getNumberActive(istrn); ++j) {
                  unsigned kind=myvals.getActiveIndex(istrn,j);
                  myvals.addDerivative( ostrn, arg_deriv_starts[k] + kind, der*myvals.getDerivative( istrn, kind ) );
              }
          }
      }
      // Clear the matrix elements for this task
      av->clearMatrixElements( myvals );
  } 
  // Turn back on matrix storing
  myvals.setMatrixStashForNormalRun(); 
  // And update the indices that there are derivatives on
  for(unsigned k=0; k<4; ++k ) {
      // This checks whether the derivatives are all on the same action or not
      bool found=false;
      for(unsigned i=0; i<k; ++i) {
          if( arg_deriv_starts[k]==arg_deriv_starts[i] ) { found=true; break; }
      }
      if( found ) continue ;
  
      unsigned istrn = getPntrToArgument(k)->getPositionInMatrixStash();
      std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( istrn ) );
      for(unsigned j=0; j<myvals.getNumberOfMatrixIndices(istrn); ++j) {
          unsigned index = arg_deriv_starts[k] + mat_indices[j];
          for(unsigned n=0; n<functions.size(); ++n) {
              unsigned ostrn = getPntrToOutput(n)->getPositionInStream(); 
              myvals.updateIndex( ostrn, index ); 
          }
      }
  } 
}

}
}
