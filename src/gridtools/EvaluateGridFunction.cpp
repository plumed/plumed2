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
#include "EvaluateGridFunction.h" 
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace gridtools {

void EvaluateGridFunction::registerKeywords( Keywords& keys ){
  keys.add("compulsory","INTERPOLATION_TYPE","spline","the method to use for interpolation.  Can be spline or floor.");
  keys.addFlag("ZERO_OUTSIDE_GRID_RANGE",false,"if we are asked to evaluate the function for a number that is outside the range of the grid set it to zero");
}

void EvaluateGridFunction::setupGridObject( Value* values, GridCoordinatesObject& gridobject ) { 
  ActionWithValue* act=values->getPntrToAction();
  std::vector<unsigned> ind( values->getRank() ), nbin( values->getRank() ); std::string gtype;
  std::vector<double> spacing( values->getRank() ), xx( values->getRank() ); std::vector<bool> pbc( values->getRank() );
  std::vector<std::string> argn( values->getRank() ), min( values->getRank() ), max( values->getRank() );
  act->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, false );
  gridobject.setup( "flat", pbc, 0, 0.0 );
}

void EvaluateGridFunction::setupGridBounds( Value* values, GridCoordinatesObject& gridobject ) {
  unsigned dimension = values->getRank();
  std::vector<std::string> argn( dimension ), min( dimension ), max( dimension ); std::string gtype;
  std::vector<unsigned> nbin( dimension ); std::vector<double> spacing( dimension ); std::vector<bool> ipbc( dimension );
  values->getPntrToAction()->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, ipbc, false );
  gridobject.setBounds( min, max, nbin, spacing );
}

void EvaluateGridFunction::read( ActionWithArguments* action ) {
  if( !action->getPntrToArgument(0)->isTimeSeries() ) {
      if( action->getPntrToArgument(0)->getRank()==0 || !action->getPntrToArgument(0)->hasDerivatives() ) action->error("should have one grid as input to this action");
  }
  // Now use this information to create a gridobject
  setupGridObject( action->getPntrToArgument(0), gridobject );
  // if( gtype=="fibonacci" ) action->error("cannot interpolate on fibonacci sphere");
  parseFlag(action,"ZERO_OUTSIDE_GRID_RANGE",set_zero_outside_range);
  if( set_zero_outside_range ) action->log.printf("  function is zero outside grid range \n");
  // Get the type of interpolation that we are doing
  std::string itype; parse(action,"INTERPOLATION_TYPE",itype);
  if( itype=="spline" ) { 
      interpolation_type=spline;
      spline_interpolator=Tools::make_unique<Interpolator>( action->getPntrToArgument(0), gridobject );
  } else if( itype=="floor" ) {
    interpolation_type=floor;
    for(unsigned j=0;j<action->getPntrToArgument(0)->getRank();++j) {
       if( !action->getPntrToArgument(0)->isTimeSeries() && gridobject.isPeriodic(j) ) action->error("floor interepolation doesn't work with periodic variables");
    }
  } else action->error("type " + itype + " of interpolation is not defined");
  action->log.printf("  generating off grid points using %s interpolation \n", itype.c_str() );
}

void EvaluateGridFunction::setup( const ActionWithArguments* action ) {
  setupGridBounds( action->getPntrToArgument(0), gridobject );
}

void EvaluateGridFunction::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  if( set_zero_outside_range && !gridobject.inbounds( args ) ) { vals[0]=0.0; return; }
  unsigned dimension = gridobject.getDimension(); plumed_dbg_assert( args.size()==dimension && vals.size()==1 ); 
  if( interpolation_type==spline ) {
      std::vector<double> der( dimension );
      vals[0] =  spline_interpolator->splineInterpolation( args, der );
      for(unsigned j=0; j<dimension; ++j) derivatives(0,j) = der[j];
  } else {
      Value* values=action->getPntrToArgument(0); std::vector<unsigned> indices(dimension); 
      gridobject.getIndices( args, indices ); unsigned nn = gridobject.getIndex(indices);
      if( values->isTimeSeries() && !gridobject.inbounds(args) ) nn = gridobject.getNbin(false)[0]-1;
      if( values->isTimeSeries() && nn==values->getShape()[0] ) {
          vals[0] = values->get( nn - 1 ); 
      } else {
          plumed_dbg_assert( nn<values->getNumberOfValues() );
          vals[0] = values->get( nn ); 
      }
      if( !values->isTimeSeries() ) {
          for(unsigned j=0; j<dimension; ++j) derivatives(0,j) = values->getGridDerivative( nn, j );
      }
  }
}

}
}
