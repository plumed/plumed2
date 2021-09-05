/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "FunctionTemplateBase.h"
#include "FunctionShortcut.h"
#include "FunctionOfVector.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace function {

class Sum : public FunctionTemplateBase {
private:
  double prefactor;
public:
  void registerKeywords( Keywords& keys ) override;
  void read( FunctionBase* action ) override;
  unsigned getRank() override;
  void setPrefactor( FunctionBase* action ) override;
  void setPeriodicityForOutputs( FunctionBase* action ) override;
  void calc( const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
};

typedef FunctionShortcut<Sum> SumShortcut;
PLUMED_REGISTER_ACTION(SumShortcut,"SUM")
PLUMED_REGISTER_ACTION(SumShortcut,"MEAN")
typedef FunctionOfVector<Sum> VectorSum;
PLUMED_REGISTER_ACTION(VectorSum,"SUM_VECTOR")
PLUMED_REGISTER_ACTION(VectorSum,"MEAN_VECTOR")

void Sum::registerKeywords( Keywords& keys ) {
  keys.use("PERIODIC");
}
 
void Sum::read( FunctionBase* action ) {
  if( action->getNumberOfArguments()!=1 ) action->error("should only be one argument to sum actions");
}

void Sum::setPrefactor( FunctionBase* action ) {
  if(action->getName().find("MEAN")!=std::string::npos) prefactor = 1. / (action->getPntrToArgument(0))->getNumberOfValues();
  else { plumed_assert( action->getName().find("SUM")!=std::string::npos ); prefactor = 1; }
}

unsigned Sum::getRank() {
  return 0;
}

void Sum::setPeriodicityForOutputs( FunctionBase* action ) {
  std::vector<std::string> period; action->parseVector("PERIODIC",period);
  if( period.size()==1 ) {
    if( period[0]!="NO") action->error("input to PERIODIC keyword does not make sense");
    action->setNotPeriodic(); return;
  } else if( period.size()!=2 ) action->error("input to PERIODIC keyword does not make sense");
  action->setPeriodic( period[0], period[1] );
}

void Sum::calc( const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  vals[0]=prefactor*args[0]; derivatives(0,0)=prefactor;
}

}
}
