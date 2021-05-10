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
#include "Function.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/ActionSetup.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION DIFFERENCE
/*
Use a switching function to determine how many of the input variables are less than a certain cutoff.

\par Examples

*/
//+ENDPLUMEDOC


class Difference :
  public Function
{
public:
  explicit Difference(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Difference,"DIFFERENCE")

void Difference::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys); keys.use("ARG");
}

Difference::Difference(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  if( arg_ends.size()!=3 ) error("difference can only take two arguments as input");
  if( getPntrToArgument(0)->isPeriodic() ) {
      if( !getPntrToArgument(1)->getPntrToAction() ) {
          if( !getPntrToArgument(1)->isPeriodic() ) error("period for input variables should be the same");
      } else {
          ActionSetup* as=dynamic_cast<ActionSetup*>( getPntrToArgument(1)->getPntrToAction() );
          if( !as && !getPntrToArgument(1)->isPeriodic() ) error("period for input variables should be the same"); 
          if( !as ) {
              std::string min0, max0; getPntrToArgument(0)->getDomain( min0, max0 );
              std::string min1, max1; getPntrToArgument(1)->getDomain( min1, max1 );
              if( min0!=min0 || max0!=max1 ) error("domain for input variables should be the same");
          }
      }
      getPeriodFromArg=0;
  } else if( getPntrToArgument(1)->isPeriodic() ) {
      ActionSetup* as=dynamic_cast<ActionSetup*>( getPntrToArgument(0)->getPntrToAction() );
      if( !as ) {
          error("period for input variables should be the same");
      } else {
          std::string min0, max0; getPntrToArgument(1)->getDomain( min0, max0 );
          getPntrToArgument(0)->setDomain( min0, max0 );
      }
      getPeriodFromArg=1;
  }
  addValueWithDerivatives(); checkRead();
}

void Difference::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  plumed_dbg_assert( args.size()==2 ); 
  addValue( 0, getPntrToArgument(0)->difference( args[1], args[0] ), myvals ); 
  addDerivative( 0, 0, 1.0, myvals ); addDerivative( 0, 1, -1.0, myvals );
}

}
}


