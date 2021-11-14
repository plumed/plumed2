/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "function/FunctionTemplateBase.h"
#include "FunctionOfGrid.h"

//+PLUMEDOC GRIDANALYSIS CUMULATIVE_INTEGRAL
/*
Calculate a cumulative integral of a function using the trapesium rule.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class CumulativeIntegral : public function::FunctionTemplateBase {
private:
  double spacing;
public:
  void registerKeywords( Keywords& keys ) override {}
  void read( ActionWithArguments* action ) override;
  bool doWithTasks() const override { return false; }
  void setPrefactor( ActionWithArguments* action, const double pref ) override;
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
};

typedef FunctionOfGrid<CumulativeIntegral> CumInt;
PLUMED_REGISTER_ACTION(CumInt,"CUMULATIVE_INTEGRAL")

void CumulativeIntegral::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=1 ) action->error("should only be one argument to this action");
}

void CumulativeIntegral::setPrefactor( ActionWithArguments* action, const double pref ) {
  spacing=pref;
} 

void CumulativeIntegral::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  plumed_assert( args.size()==vals.size() ); vals[0]=0;
  for(unsigned i=1; i<vals.size(); ++i) vals[i] = vals[i-1] + (args[i-1]+args[i])*spacing / 2.;
}

}
}
