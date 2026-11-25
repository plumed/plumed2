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
#include "core/ActionWithVector.h"
#include "core/ActionRegister.h"
#include "tools/LeptonCall.h"
#include "tools/Angle.h"

namespace PLMD {
namespace symfunc {

//+PLUMEDOC COLVAR GSYMFUNC_THREEBODY
/*
Calculate functions of the coordinates of the coordinates of all pairs of bonds in the first coordination sphere of an atom

\par Examples

*/
//+ENDPLUMEDOC

class ThreeBodyGFunctions : public ActionWithVector {
private:
  std::vector<LeptonCall> functions;
public:
  static void registerKeywords( Keywords& keys );
  explicit ThreeBodyGFunctions(const ActionOptions&);
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  void calculate() override ;
  unsigned getNumberOfDerivatives() override;
  void performTask( const unsigned& task_index, MultiValue& myvals ) const override ;
};

PLUMED_REGISTER_ACTION(ThreeBodyGFunctions,"GSYMFUNC_THREEBODY")

void ThreeBodyGFunctions::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.use("ARG");
  keys.add("compulsory","WEIGHT","the matrix that contains the weights that should be used for each connection");
  keys.add("numbered","FUNCTION","the parameters of the function you would like to compute");
  ActionWithValue::useCustomisableComponents( keys );
}

ThreeBodyGFunctions::ThreeBodyGFunctions(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao) {
  if( getNumberOfArguments()!=3 ) {
    error("found wrong number of arguments in input");
  }
  std::vector<Value*> wval;
  parseArgumentList("WEIGHT",wval);
  if( wval.size()!=1 ) {
    error("keyword WEIGHT should be provided with the label of a single action");
  }

  for(unsigned i=0; i<3; ++i) {
    if( getPntrToArgument(i)->getRank()!=2 ) {
      error("input argument should be a matrix");
    }
    if( wval[0]->getShape()[0]!=getPntrToArgument(i)->getShape()[0] || wval[0]->getShape()[1]!=getPntrToArgument(i)->getShape()[1] ) {
      error("mismatched shapes of matrices in input");
    }
  }
  log.printf("  using bond weights from matrix labelled %s \n",wval[0]->getName().c_str() );
  // Rerequest the arguments
  std::vector<Value*> myargs( getArguments() );
  myargs.push_back( wval[0] );
  requestArguments( myargs );
  for(unsigned i=0; i<myargs.size(); ++i) {
    myargs[i]->buildDataStore();
  }
  std::vector<unsigned> shape(1);
  shape[0] = getPntrToArgument(0)->getShape()[0];

  // And now read the functions to compute
  for(int i=1;; i++) {
    std::string myfunc, mystr, lab, num;
    Tools::convert(i,num);
    if( !parseNumbered("FUNCTION",i,mystr ) ) {
      break;
    }
    std::vector<std::string> data=Tools::getWords(mystr);
    if( !Tools::parse(data,"LABEL",lab ) ) {
      error("found no LABEL in FUNCTION" + num + " specification");
    }
    addComponent( lab, shape );
    componentIsNotPeriodic( lab );
    if( !Tools::parse(data,"FUNC",myfunc) ) {
      error("found no FUNC in FUNCTION" + num + " specification");
    }
    log.printf("  component labelled %s is computed using %s \n",lab.c_str(), myfunc.c_str() );
    functions.push_back( LeptonCall() );
    std::vector<std::string> argnames(1);
    argnames[0]="ajik";
    if( myfunc.find("rij")!=std::string::npos ) {
      argnames.push_back("rij");
    }
    if( myfunc.find("rik")!=std::string::npos ) {
      if( argnames.size()<2 ) {
        error("if you have a function of rik it must also be a function of rij -- email gareth.tribello@gmail.com if this is a problem");
      }
      argnames.push_back("rik");
    }
    if( myfunc.find("rjk")!=std::string::npos ) {
      if( argnames.size()<2 ) {
        error("if you have a function of rjk it must also be a function of rij and rik -- email gareth.tribello@gmail.com if this is a problem");
      }
      argnames.push_back("rjk");
    }
    functions[i-1].set( myfunc, argnames, this, true );
  }
  checkRead();
}

std::string ThreeBodyGFunctions::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    if( getConstPntrToComponent(i)->getName().find(cname)!=std::string::npos ) {
      std::string num;
      Tools::convert( i+1, num );
      return "the function defined by the FUNCTION" + num + " keyword";
    }
  }
  plumed_error();
  return "";
}

unsigned ThreeBodyGFunctions::getNumberOfDerivatives() {
  return 0;
}

void ThreeBodyGFunctions::calculate() {
  runAllTasks();
}

void ThreeBodyGFunctions::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  const Value* wval = getPntrToArgument(3);
  const Value* xval = getPntrToArgument(0);
  const Value* yval = getPntrToArgument(1);
  const Value* zval = getPntrToArgument(2);
  Angle angle;
  Vector disti, distj;
  unsigned matsize = wval->getNumberOfValues();
  std::vector<double> values(4);
  std::vector<Vector> der_i(4), der_j(4);
  unsigned nbonds = wval->getRowLength( task_index ), ncols = wval->getShape()[1];
  for(unsigned i=0; i<nbonds; ++i) {
    unsigned ipos = ncols*task_index + wval->getRowIndex( task_index, i );
    double weighti = wval->get( ipos );
    if( weighti<epsilon ) {
      continue ;
    }
    disti[0] = xval->get( ipos );
    disti[1] = yval->get( ipos );
    disti[2] = zval->get( ipos );
    values[1] = disti.modulo2();
    der_i[1]=2*disti;
    der_i[2].zero();
    for(unsigned j=0; j<i; ++j) {
      unsigned jpos = ncols*task_index + wval->getRowIndex( task_index, j );
      double weightj = wval->get( jpos );
      if( weightj<epsilon ) {
        continue ;
      }
      distj[0] = xval->get( jpos );
      distj[1] = yval->get( jpos );
      distj[2] = zval->get( jpos );
      values[2] = distj.modulo2();
      der_j[1].zero();
      der_j[2]=2*distj;
      der_i[3] = ( disti - distj );
      values[3] = der_i[3].modulo2();
      der_i[3] = 2*der_i[3];
      der_j[3] = -der_i[3];
      // Compute angle between bonds
      values[0] = angle.compute( disti, distj, der_i[0], der_j[0] );
      // Compute product of weights
      double weightij = weighti*weightj;
      // Now compute all symmetry functions
      for(unsigned n=0; n<functions.size(); ++n) {
        unsigned ostrn = getConstPntrToComponent(n)->getPositionInStream();
        double nonweight = functions[n].evaluate( values );
        myvals.addValue( ostrn, nonweight*weightij );
        if( doNotCalculateDerivatives() ) {
          continue;
        }

        for(unsigned m=0; m<functions[n].getNumberOfArguments(); ++m) {
          double der = weightij*functions[n].evaluateDeriv( m, values );
          myvals.addDerivative( ostrn, ipos, der*der_i[m][0] );
          myvals.addDerivative( ostrn, matsize+ipos, der*der_i[m][1] );
          myvals.addDerivative( ostrn, 2*matsize+ipos, der*der_i[m][2] );
          myvals.addDerivative( ostrn, jpos, der*der_j[m][0] );
          myvals.addDerivative( ostrn, matsize+jpos, der*der_j[m][1] );
          myvals.addDerivative( ostrn, 2*matsize+jpos, der*der_j[m][2] );
        }
        myvals.addDerivative( ostrn, 3*matsize+ipos, nonweight*weightj );
        myvals.addDerivative( ostrn, 3*matsize+jpos, nonweight*weighti );
      }
    }
  }
  if( doNotCalculateDerivatives() ) {
    return ;
  }

  // And update the elements that have derivatives
  // Needs a separate loop here as there may be forces from j
  for(unsigned i=0; i<nbonds; ++i) {
    unsigned ipos = ncols*task_index + wval->getRowIndex( task_index, i );
    double weighti = wval->get( ipos );
    if( weighti<epsilon ) {
      continue ;
    }

    for(unsigned n=0; n<functions.size(); ++n) {
      unsigned ostrn = getConstPntrToComponent(n)->getPositionInStream();
      myvals.updateIndex( ostrn, ipos );
      myvals.updateIndex( ostrn, matsize+ipos );
      myvals.updateIndex( ostrn, 2*matsize+ipos );
      myvals.updateIndex( ostrn, 3*matsize+ipos );
    }
  }
}

}
}
