/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "MoreThan.h"
#include "FunctionShortcut.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"

#include <cmath>

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION MORE_THAN
/*
Use a switching function to determine how many of the input variables are more than a certain cutoff.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION MORE_THAN_VECTOR
/*
Use a switching function to determine how many of elements in the input vector are more than a certain cutoff.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR MORE_THAN_MATRIX
/*
Transform all the elements of a matrix using a switching function that is one when the input value is larger than a threshold

\par Examples

*/
//+ENDPLUMEDOC

typedef FunctionShortcut<MoreThan> MoreThanShortcut;
PLUMED_REGISTER_ACTION(MoreThanShortcut,"MORE_THAN")
typedef FunctionOfVector<MoreThan> VectorMoreThan;
PLUMED_REGISTER_ACTION(VectorMoreThan,"MORE_THAN_VECTOR")
typedef FunctionOfMatrix<MoreThan> MatrixMoreThan;
PLUMED_REGISTER_ACTION(MatrixMoreThan,"MORE_THAN_MATRIX")

void MoreThan::registerKeywords(Keywords& keys) {
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.addFlag("SQUARED",false,"is the input quantity the square of the value that you would like to apply the switching function to");
  keys.setValueDescription("scalar/vector/matrix","a function that is one if the if the input is more than a threshold");
}

void MoreThan::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=1 ) {
    action->error("should only be one argument to more_than actions");
  }
  if( action->getPntrToArgument(0)->isPeriodic() ) {
    action->error("cannot use this function on periodic functions");
  }


  std::string sw,errors;
  action->parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
    if( errors.length()!=0 ) {
      action->error("problem reading SWITCH keyword : " + errors );
    }
  } else {
    int nn=6;
    int mm=0;
    double d0=0.0;
    double r0=0.0;
    action->parse("R_0",r0);
    if(r0<=0.0) {
      action->error("R_0 should be explicitly specified and positive");
    }
    action->parse("D_0",d0);
    action->parse("NN",nn);
    action->parse("MM",mm);
    switchingFunction.set(nn,mm,r0,d0);
  }
  action->log<<"  using switching function with cutoff "<<switchingFunction.description()<<"\n";
  action->parseFlag("SQUARED",squared);
  if( squared ) {
    action->log<<"  input quantity is square of quantity that switching function acts upon\n";
  }
}

void MoreThan::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  plumed_dbg_assert( args.size()==1 );
  if( squared ) {
    vals[0] = 1.0 - switchingFunction.calculateSqr( args[0], derivatives(0,0) );
  } else {
    vals[0] = 1.0 - switchingFunction.calculate( args[0], derivatives(0,0) );
  }
  derivatives(0,0) = -args[0]*derivatives(0,0);
}

}
}


