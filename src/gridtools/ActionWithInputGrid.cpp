/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "ActionWithInputGrid.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace gridtools {

void ActionWithInputGrid::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  EvaluateGridFunction ii; ii.registerKeywords( keys );  
}

ActionWithInputGrid::ActionWithInputGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  firststep(true),
  set_zero_outside_range(false)
{
  if( getNumberOfArguments()!=1 ) error("should be exactly one argument to this action");
  if( getPntrToArgument(0)->getRank()==0 || (!getPntrToArgument(0)->isTimeSeries() && !getPntrToArgument(0)->hasDerivatives()) ) error("input to this action should be a grid");

  my_interpolator.read( this );
}

double ActionWithInputGrid::getFunctionValueAndDerivatives( const std::vector<double>& x, std::vector<double>& der ) const {
  std::vector<double> value(1); Matrix<double> derivatives( 1, der.size() );
  my_interpolator.calc( this, x, value, derivatives );
  for(unsigned i=0; i<der.size(); ++i) der[i] = derivatives(0,i);
  return value[0]; 
}

void ActionWithInputGrid::setupGridObject() {
  plumed_assert( firststep ); 
  my_interpolator.setup( this );
}

void ActionWithInputGrid::doTheCalculation() {
  if( getPntrToArgument(0)->getShape()[0]==0 ) return;
  if( firststep ) { setupGridObject(); firststep=false; finishOutputSetup(); }
  if( getFullNumberOfTasks()>0 ) { runAllTasks(); jobsAfterLoop(); }
  else runTheCalculation();
}

void ActionWithInputGrid::calculate() {
  if( actionInChain() ) return ;
  doTheCalculation();
}

void ActionWithInputGrid::update() {
  if( skipUpdate() ) return;
  doTheCalculation();
}

void ActionWithInputGrid::runFinalJobs() {
  if( skipUpdate() ) return;
  doTheCalculation();
}

}
}

