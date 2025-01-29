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
#include "function/FunctionShortcut.h"
#include "function/FunctionOfScalar.h"
#include "function/FunctionOfVector.h"
#include "core/ActionRegister.h"
#include "function/FunctionTemplateBase.h"

#include <cmath>

namespace PLMD {
namespace refdist {

//+PLUMEDOC FUNCTION DIFFERENCE
/*
Calculate the differences between two scalars

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION DIFFERENCE_SCALAR
/*
Calculate the differences between two scalars

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION DIFFERENCE_VECTOR
/*
Calculate the differences between the elements of two vectors

\par Examples

*/
//+ENDPLUMEDOC


class Difference : public function::FunctionTemplateBase {
private:
  bool periodic;
  std::string min0, max0;
public:
  void registerKeywords(Keywords& keys) override ;
  void read( ActionWithArguments* action ) override;
  void setPeriodicityForOutputs( ActionWithValue* action ) override;
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
};


typedef function::FunctionShortcut<Difference> DifferenceShortcut;
PLUMED_REGISTER_ACTION(DifferenceShortcut,"DIFFERENCE")
typedef function::FunctionOfScalar<Difference> ScalarDifference;
PLUMED_REGISTER_ACTION(ScalarDifference,"DIFFERENCE_SCALAR")
typedef function::FunctionOfVector<Difference> VectorDifference;
PLUMED_REGISTER_ACTION(VectorDifference,"DIFFERENCE_VECTOR")

void Difference::registerKeywords(Keywords& keys) {
  keys.setValueDescription("scalar/vector","a function that measures the difference");
}

void Difference::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=2 ) {
    action->error("should be two arguments to this action");
  }
  if( action->getPntrToArgument(0)->getRank()==action->getPntrToArgument(1)->getRank() ) {
    std::vector<unsigned> shape( action->getPntrToArgument(0)->getShape() );
    for(unsigned i=0; i<shape.size(); ++i) {
      if( shape[i]!=action->getPntrToArgument(1)->getShape()[i] ) {
        action->error("shapes of input actions do not match");
      }
    }
  }

  periodic=false;
  if( action->getPntrToArgument(0)->isPeriodic() ) {
    periodic=true;
    action->getPntrToArgument(0)->getDomain( min0, max0 );
    if( !action->getPntrToArgument(1)->isConstant() && !action->getPntrToArgument(1)->isPeriodic() ) {
      action->error("period for input variables " + action->getPntrToArgument(0)->getName() + " and " + action->getPntrToArgument(1)->getName() + " should be the same 0");
    }
    if( !action->getPntrToArgument(1)->isConstant() ) {
      std::string min1, max1;
      action->getPntrToArgument(1)->getDomain( min1, max1 );
      if( min0!=min0 || max0!=max1 ) {
        action->error("domain for input variables should be the same");
      }
    } else {
      action->getPntrToArgument(1)->setDomain( min0, max0 );
    }
  } else if( action->getPntrToArgument(1)->isPeriodic() ) {
    periodic=true;
    action->getPntrToArgument(1)->getDomain( min0, max0 );
    if( !action->getPntrToArgument(1)->isConstant() ) {
      action->error("period for input variables " + action->getPntrToArgument(0)->getName() + " and " + action->getPntrToArgument(1)->getName() + " should be the same 1");
    } else {
      action->getPntrToArgument(0)->setDomain( min0, max0 );
    }
  }
}

void Difference::setPeriodicityForOutputs( ActionWithValue* action ) {
  if( periodic ) {
    action->setPeriodic( min0, max0 );
  } else {
    action->setNotPeriodic();
  }
}

void Difference::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  plumed_dbg_assert( args.size()==2 );
  vals[0] = action->getPntrToArgument(0)->difference( args[1], args[0] );
  derivatives(0,0) = 1.0;
  derivatives(0,1)=-1;
}

}
}


