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

namespace PLMD {
namespace gridtools {

void ActionWithInputGrid::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","INTERPOLATION_TYPE","spline","the method to use for interpolation.  Can be spline or floor.");
}

ActionWithInputGrid::ActionWithInputGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  my_interpolator( getPntrToArgument(0), gridobject ),
  firststep(true),
  set_zero_outside_range(false)
{
  if( getNumberOfArguments()!=1 ) error("should be exactly one argument to this action");
  if( getPntrToArgument(0)->getRank()==0 || (!getPntrToArgument(0)->isTimeSeries() && !getPntrToArgument(0)->hasDerivatives()) ) error("input to this action should be a grid");

  unsigned dimension = getPntrToArgument(0)->getRank();
  std::vector<std::string> argn( dimension ), min( dimension ), max( dimension ); std::string gtype;
  std::vector<unsigned> nbin( dimension ); std::vector<double> spacing( dimension ); std::vector<bool> ipbc( dimension );
  (getPntrToArgument(0)->getPntrToAction())->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, ipbc, false );
  if( gtype=="flat" ) gridobject.setup( "flat", ipbc, 0, 0.0 );
  else if( gtype=="fibonacci" ) gridobject.setup( "fibonacci", ipbc, nbin[0], spacing[0] );
  else plumed_merror("unknown grid type");

  std::string itype; parse("INTERPOLATION_TYPE",itype);
  if( itype=="spline" ) interpolation_type=spline;
  else if( itype=="floor" ) {
    interpolation_type=floor;
    for(unsigned j=0;j<dimension;++j) if( !getPntrToArgument(0)->isTimeSeries() && ipbc[j] ) error("floor interepolation doesn't work with periodic variables");
  } else error("type " + itype + " of interpolation is not defined");
  log.printf("  generating off grid points using %s interpolation \n", itype.c_str() );
}

double ActionWithInputGrid::getFunctionValueAndDerivatives( const std::vector<double>& x, std::vector<double>& der ) const {
  plumed_dbg_assert( gridobject.getGridType()=="flat" ); 
  if( set_zero_outside_range && !gridobject.inbounds( x ) ) return 0.0;
  double value;

  // loop over neighbors
  if( interpolation_type==spline ) {
      value = my_interpolator.splineInterpolation( x, der );
  } else if ( interpolation_type==floor ) {
      unsigned dimension = gridobject.getDimension(); std::vector<unsigned> indices(dimension); 
      gridobject.getIndices( x, indices ); unsigned nn = gridobject.getIndex(indices); 
      if( getPntrToArgument(0)->isTimeSeries() && !gridobject.inbounds(x) ) nn = gridobject.getNbin(false)[0]-1;
      value = getFunctionValue( nn ); if( getPntrToArgument(0)->isTimeSeries() ) return value;
      for(unsigned j=0; j<dimension; ++j) der[j] = getPntrToArgument(0)->getGridDerivative( nn, j ); 
  }
  return value;
}

void ActionWithInputGrid::setupGridObject() {
  plumed_assert( firststep ); 
  if( gridobject.getGridType()=="flat" ) {
    unsigned dimension = getPntrToArgument(0)->getRank();
    std::vector<std::string> argn( dimension ), min( dimension ), max( dimension ); std::string gtype;
    std::vector<unsigned> nbin( dimension ); std::vector<double> spacing( dimension ); std::vector<bool> ipbc( dimension );
    (getPntrToArgument(0)->getPntrToAction())->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, ipbc, false );
    gridobject.setBounds( min, max, nbin, spacing );
  }
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

