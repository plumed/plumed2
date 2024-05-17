/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionTemplateBase.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION HIGHEST
/*
This function can be used to find the highest colvar by magnitude in a set.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR HIGHEST_SCALAR
/*
Calculate the highest of a set of sclalar arguments

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR HIGHEST_VECTOR
/*
Calculate the largest element in a vector of inputs

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION LOWEST
/*
This function can be used to find the lowest colvar by magnitude in a set.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR LOWEST_SCALAR
/*
Calculate the lowest of a set of sclalar arguments

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR LOWEST_VECTOR
/*
Calculate the lowest element in a vector of inputs

\par Examples

*/
//+ENDPLUMEDOC

class Highest : public FunctionTemplateBase {
private:
  bool min, scalar_out;
public:
  void registerKeywords( Keywords& keys ) override ;
  void read( ActionWithArguments* action ) override;
  bool zeroRank() const override { return scalar_out; }
  bool doWithTasks() const override { return !scalar_out; }
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
};

typedef FunctionShortcut<Highest> HighestShortcut;
PLUMED_REGISTER_ACTION(HighestShortcut,"HIGHEST")
PLUMED_REGISTER_ACTION(HighestShortcut,"LOWEST")
typedef FunctionOfScalar<Highest> ScalarHighest;
PLUMED_REGISTER_ACTION(ScalarHighest,"HIGHEST_SCALAR")
PLUMED_REGISTER_ACTION(ScalarHighest,"LOWEST_SCALAR")
typedef FunctionOfVector<Highest> VectorHighest;
PLUMED_REGISTER_ACTION(VectorHighest,"HIGHEST_VECTOR")
PLUMED_REGISTER_ACTION(VectorHighest,"LOWEST_VECTOR")

void Highest::registerKeywords( Keywords& keys ) {
  if( keys.getDisplayName().find("LOWEST") ) keys.setValueDescription("the lowest of the input values");
  else keys.setValueDescription("the highest of the input values");
}

void Highest::read( ActionWithArguments* action ) {
  min=action->getName().find("LOWEST")!=std::string::npos; if( !min ) plumed_assert( action->getName().find("HIGHEST")!=std::string::npos );
  for(unsigned i=0; i<action->getNumberOfArguments(); ++i) {
    if( action->getPntrToArgument(i)->isPeriodic() ) action->error("Cannot sort periodic values (check argument "+ action->getPntrToArgument(i)->getName() +")");
  }
  scalar_out = action->getNumberOfArguments()==1;
  if( scalar_out && action->getPntrToArgument(0)->getRank()==0 ) action->error("sorting a single scalar is trivial");
}

void Highest::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  if( min ) {
    vals[0] = *std::min_element(args.begin(), args.end());
    derivatives(0,std::min_element(args.begin(), args.end()) - args.begin()) = 1;
  } else {
    vals[0] = *std::max_element(args.begin(), args.end());
    derivatives(0,std::max_element(args.begin(), args.end()) - args.begin()) = 1;
  }
}

}
}


