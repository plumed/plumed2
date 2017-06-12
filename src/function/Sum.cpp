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
#include "tools/SwitchingFunction.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION SUM
/*
Calculate the sum of the input values

\par Examples

*/
//+ENDPLUMEDOC


class Sum :
  public Function
{
public:
  explicit Sum(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Sum,"SUM")

void Sum::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG"); 
}

Sum::Sum(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  for(unsigned i=0;i<getNumberOfComponents();++i){
     if( getPntrToComponent(i)->isPeriodic() ) error("cannot use this function on periodic functions");
  }

  std::vector<unsigned> shape(0);
  if( arg_ends[1]-arg_ends[0]==1 ){
      ActionWithValue::addValueWithDerivatives( shape ); setNotPeriodic();
  } else {
      std::string num;
      for(unsigned i=0;i<arg_ends.size()-1;++i){ 
         Tools::convert(i+1,num); ActionWithValue::addComponentWithDerivatives( "arg_" + num, shape ); 
         componentIsNotPeriodic( "arg_" + num );
      }
  }
  checkRead();
}

void Sum::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  plumed_dbg_assert( args.size()==1 ); setValue( 0, args[0], myvals ); addDerivative( 0, 0, 1.0, myvals );
}

}
}


