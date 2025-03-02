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
#include "core/ParallelTaskManager.h"
#include "core/MatrixView.h"
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

class ThreeBodyGFunctionsInput {
public:
  bool multi_action_input;
  std::vector<std::string> funcstr;
  std::vector<LeptonCall> functions;
  std::vector<MatrixView> args;
  ThreeBodyGFunctionsInput& operator=( const ThreeBodyGFunctionsInput& m ) {
    args.resize(4);
    multi_action_input = m.multi_action_input;
    for(unsigned i=0; i<m.funcstr.size(); ++i) addFunction( m.funcstr[i], NULL ); 
    return *this;
  }
  void addFunction( std::string myfunc, ActionWithVector* action ) {
    funcstr.push_back( myfunc );
    functions.push_back( LeptonCall() );
    std::vector<std::string> argnames(1);
    argnames[0]="ajik";
    if( myfunc.find("rij")!=std::string::npos ) {
      argnames.push_back("rij");
    }   
    if( myfunc.find("rik")!=std::string::npos ) {
      if( action && argnames.size()<2 ) { 
        action->error("if you have a function of rik it must also be a function of rij -- email gareth.tribello@gmail.com if this is a problem");
      }
      argnames.push_back("rik");
    }
    if( myfunc.find("rjk")!=std::string::npos ) {
      if( action && argnames.size()<2 ) {
        action->error("if you have a function of rjk it must also be a function of rij and rik -- email gareth.tribello@gmail.com if this is a problem");
      }
      argnames.push_back("rjk");
    }
    functions[functions.size()-1].set( myfunc, argnames, action, true );
  }
};

class ThreeBodyGFunctions : public ActionWithVector {
private:
  ParallelTaskManager<ThreeBodyGFunctions,ThreeBodyGFunctionsInput> taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  explicit ThreeBodyGFunctions(const ActionOptions&);
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  void calculate() override ;
  unsigned getNumberOfDerivatives() override;
  void performTask( const unsigned& task_index, MultiValue& myvals ) const override ;
  static void performTask( std::size_t task_index, const ThreeBodyGFunctionsInput& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output );
};

PLUMED_REGISTER_ACTION(ThreeBodyGFunctions,"GSYMFUNC_THREEBODY")

void ThreeBodyGFunctions::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","matrix","three matrices containing the bond vectors of interest");
  keys.addInputKeyword("compulsory","WEIGHT","matrix","the matrix that contains the weights that should be used for each connection");
  keys.add("numbered","FUNCTION","the parameters of the function you would like to compute");
  ParallelTaskManager<ThreeBodyGFunctions,ThreeBodyGFunctionsInput>::registerKeywords( keys );
  ActionWithValue::useCustomisableComponents( keys );
}

ThreeBodyGFunctions::ThreeBodyGFunctions(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this) {
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
  std::vector<std::size_t> shape(1);
  shape[0] = getPntrToArgument(0)->getShape()[0];

  // And now read the functions to compute
  ThreeBodyGFunctionsInput input;
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
    input.addFunction( myfunc, this ); 
  }
  checkRead();
  input.multi_action_input = getPntrToArgument(3)->getPntrToAction()!=getPntrToArgument(0)->getPntrToAction();
  taskmanager.setupParallelTaskManager( 1, 0, 0 );
  taskmanager.setActionInput( input );
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
  std::size_t s = 0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    taskmanager.getActionInput().args[i].setup( s, getPntrToArgument(i) );
    s += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  taskmanager.runAllTasks();
}

void ThreeBodyGFunctions::performTask( std::size_t task_index, const ThreeBodyGFunctionsInput& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output ) {
  const double* xpntr=NULL; 
  const double* ypntr=NULL; 
  const double* zpntr=NULL;
  std::vector<double> xvals, yvals, zvals;
  std::size_t rowlen = actiondata.args[3].bookeeping[(1+actiondata.args[3].ncols)*task_index];
  View<const double,helpers::dynamic_extent> wval( input.inputdata.data() + actiondata.args[3].start + actiondata.args[3].ncols*task_index, rowlen ); 
  if( actiondata.multi_action_input ) {
      xvals.resize( rowlen );
      yvals.resize( rowlen );
      zvals.resize( rowlen );
      View<const std::size_t,helpers::dynamic_extent> wbooks( actiondata.args[3].bookeeping.data()+(1+actiondata.args[3].ncols)*task_index+1, rowlen);
      for(unsigned i=0; i<rowlen; ++i) {
          xvals[i] = MatrixView::getElement( task_index, wbooks[i], actiondata.args[0], input.inputdata );
          yvals[i] = MatrixView::getElement( task_index, wbooks[i], actiondata.args[1], input.inputdata );
          zvals[i] = MatrixView::getElement( task_index, wbooks[i], actiondata.args[2], input.inputdata );
      }
      xpntr = xvals.data();
      ypntr = yvals.data();
      zpntr = zvals.data(); 
  } else {
     xpntr = input.inputdata.data() + actiondata.args[0].start + actiondata.args[0].ncols*task_index;
     ypntr = input.inputdata.data() + actiondata.args[1].start + actiondata.args[1].ncols*task_index;
     zpntr = input.inputdata.data() + actiondata.args[2].start + actiondata.args[2].ncols*task_index; 
  }
  View<const double,helpers::dynamic_extent> xval( xpntr, rowlen );
  View<const double,helpers::dynamic_extent> yval( ypntr, rowlen );
  View<const double,helpers::dynamic_extent> zval( zpntr, rowlen );

  Angle angle;
  Vector disti, distj;
  std::vector<double> values(4);
  std::vector<Vector> der_i(4), der_j(4); 
  for(unsigned i=0; i<rowlen; ++i) {
      if( wval[i]<epsilon ) {
          continue;
      }
      disti[0] = xval[i]; 
      disti[1] = yval[i]; 
      disti[2] = zval[i];  
      values[1] = disti.modulo2();
      der_i[1]=2*disti;
      der_i[2].zero();
      for(unsigned j=0; j<i; ++j) {
          if( wval[j]<epsilon) {
              continue;
          }
          distj[0] = xval[j]; 
          distj[1] = yval[j]; 
          distj[2] = zval[j]; 
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
          double weightij = wval[i]*wval[j];
          for(unsigned n=0; n<actiondata.functions.size(); ++n) {
            double nonweight = actiondata.functions[n].evaluate( values );
            output.values[n] += nonweight*weightij;

            if( input.noderiv ) {
                continue;
            }
          }
      }
  } 
}

void ThreeBodyGFunctions::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  const Value* wval = getPntrToArgument(3);
  const Value* xval = getPntrToArgument(0);
  const Value* yval = getPntrToArgument(1);
  const Value* zval = getPntrToArgument(2);
  Angle angle;
  Vector disti, distj;
  unsigned matsize = wval->getNumberOfStoredValues();
  std::vector<double> values(4);
  std::vector<Vector> der_i(4), der_j(4);
  unsigned nbonds = wval->getRowLength( task_index ), ncols = wval->getNumberOfColumns();
  if( taskmanager.getActionInput().multi_action_input ) {
    matsize = wval->getNumberOfValues();
    ncols = wval->getShape()[1];
  }
  for(unsigned i=0; i<nbonds; ++i) {
    unsigned ipos = ncols*task_index + i;  //wval->getRowIndex( task_index, i );
    if( taskmanager.getActionInput().multi_action_input ) {
      ipos = ncols*task_index + wval->getRowIndex( task_index, i );
    }
    double weighti = wval->get( ipos, taskmanager.getActionInput().multi_action_input );
    if( weighti<epsilon ) {
      continue ;
    }
    disti[0] = xval->get( ipos, taskmanager.getActionInput().multi_action_input );
    disti[1] = yval->get( ipos, taskmanager.getActionInput().multi_action_input );
    disti[2] = zval->get( ipos, taskmanager.getActionInput().multi_action_input );
    values[1] = disti.modulo2();
    der_i[1]=2*disti;
    der_i[2].zero();
    for(unsigned j=0; j<i; ++j) {
      unsigned jpos = ncols*task_index + j;  // wval->getRowIndex( task_index, j );
      if( taskmanager.getActionInput().multi_action_input ) {
        jpos = ncols*task_index + wval->getRowIndex( task_index, j );
      }
      double weightj = wval->get( jpos, taskmanager.getActionInput().multi_action_input );
      if( weightj<epsilon ) {
        continue ;
      }
      distj[0] = xval->get( jpos, taskmanager.getActionInput().multi_action_input );
      distj[1] = yval->get( jpos, taskmanager.getActionInput().multi_action_input );
      distj[2] = zval->get( jpos, taskmanager.getActionInput().multi_action_input );
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
      for(unsigned n=0; n<taskmanager.getActionInput().functions.size(); ++n) {
        double nonweight = taskmanager.getActionInput().functions[n].evaluate( values );
        myvals.addValue( n, nonweight*weightij );
        if( doNotCalculateDerivatives() ) {
          continue;
        }

        for(unsigned m=0; m<taskmanager.getActionInput().functions[n].getNumberOfArguments(); ++m) {
          double der = weightij*taskmanager.getActionInput().functions[n].evaluateDeriv( m, values );
          myvals.addDerivative( n, ipos, der*der_i[m][0] );
          myvals.addDerivative( n, matsize+ipos, der*der_i[m][1] );
          myvals.addDerivative( n, 2*matsize+ipos, der*der_i[m][2] );
          myvals.addDerivative( n, jpos, der*der_j[m][0] );
          myvals.addDerivative( n, matsize+jpos, der*der_j[m][1] );
          myvals.addDerivative( n, 2*matsize+jpos, der*der_j[m][2] );
        }
        myvals.addDerivative( n, 3*matsize+ipos, nonweight*weightj );
        myvals.addDerivative( n, 3*matsize+jpos, nonweight*weighti );
      }
    }
  }
  if( doNotCalculateDerivatives() ) {
    return ;
  }

  // And update the elements that have derivatives
  // Needs a separate loop here as there may be forces from j
  for(unsigned i=0; i<nbonds; ++i) {
    unsigned ipos = ncols*task_index + i; // wval->getRowIndex( task_index, i );
    double weighti = wval->get( ipos, false );
    if( weighti<epsilon ) {
      continue ;
    }

    for(unsigned n=0; n<taskmanager.getActionInput().functions.size(); ++n) {
      myvals.updateIndex( n, ipos );
      myvals.updateIndex( n, matsize+ipos );
      myvals.updateIndex( n, 2*matsize+ipos );
      myvals.updateIndex( n, 3*matsize+ipos );
    }
  }
}

}
}
