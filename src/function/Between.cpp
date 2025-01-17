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
#include "Between.h"
#include "FunctionShortcut.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"

#include <cmath>

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION BETWEEN
/*
Use a switching function to determine how many of the input variables are within a certain range.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION BETWEEN_VECTOR
/*
Use a switching function to determine how many of the input components are within a certain range

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR BETWEEN_MATRIX
/*
Transform all the elements of a matrix using a switching function that is oen when the input value is within a particular range

\par Examples

*/
//+ENDPLUMEDOC

typedef FunctionShortcut<Between> BetweenShortcut;
PLUMED_REGISTER_ACTION(BetweenShortcut,"BETWEEN")
typedef FunctionOfVector<Between> VectorBetween;
PLUMED_REGISTER_ACTION(VectorBetween,"BETWEEN_VECTOR")
typedef FunctionOfMatrix<Between> MatrixBetween;
PLUMED_REGISTER_ACTION(MatrixBetween,"BETWEEN_MATRIX")

void Between::registerKeywords(Keywords& keys) {
  keys.add("compulsory","LOWER","the lower boundary for this particular bin");
  keys.add("compulsory","UPPER","the upper boundary for this particular bin");
  keys.add("compulsory","SMEAR","0.5","the ammount to smear the Gaussian for each value in the distribution");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous function defined above. "
           "The following provides information on the \\ref histogrambead that are available. "
           "When this keyword is present you no longer need the LOWER, UPPER, SMEAR and KERNEL keywords.");
  keys.setValueDescription("a function that is one if the input falls within a particular range and zero otherwise");
}

void Between::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=1 ) {
    action->error("should only be one argument to between actions");
  }

  std::string str_min, str_max, tstr_min, tstr_max;
  bool isPeriodic = action->getPntrToArgument(0)->isPeriodic();
  if( isPeriodic ) {
    action->getPntrToArgument(0)->getDomain( str_min, str_max );
  }

  std::string hinput;
  action->parse("SWITCH",hinput);
  if(hinput.length()==0) {
    std::string low, up, sme;
    action->parse("LOWER",low);
    action->parse("UPPER",up);
    action->parse("SMEAR",sme);
    hinput = "GAUSSIAN LOWER=" + low + " UPPER=" + up + " SMEAR=" + sme;
  }
  std::string errors;
  hist.set( hinput, errors );
  if( errors.size()!=0 ) {
    action->error( errors );
  }
  action->log.printf("  %s \n", hist.description().c_str() );

  if( !isPeriodic ) {
    hist.isNotPeriodic();
  } else {
    double min;
    Tools::convert( str_min, min );
    double max;
    Tools::convert( str_max, max );
    hist.isPeriodic( min, max );
  }
}

void Between::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  plumed_dbg_assert( args.size()==1 );
  vals[0] = hist.calculate( args[0], derivatives(0,0) );
}

}
}


