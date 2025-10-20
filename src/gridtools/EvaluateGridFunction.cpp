/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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

void EvaluateGridFunction::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","INTERPOLATION_TYPE","spline","the method to use for interpolation.  Can be spline, linear, ceiling or floor.");
  keys.addFlag("ZERO_OUTSIDE_GRID_RANGE",false,"if we are asked to evaluate the function for a number that is outside the range of the grid set it to zero");
  keys.setValueDescription("scalar","interpolation of the input grid to get the value of the function at the input arguments");
}

std::vector<bool> EvaluateGridFunction::getPbc() const {
  const GridCoordinatesObject & gridobject = getGridObject();
  std::vector<bool> ipbc( gridobject.getDimension() );
  for(unsigned i=0; i<ipbc.size(); ++i) {
    ipbc[i] = gridobject.isPeriodic(i);
  }
  return ipbc;
}

void EvaluateGridFunction::read( EvaluateGridFunction& func, ActionWithArguments* action, function::FunctionOptions& options ) {
  func.function = action->getPntrToArgument(0);
  if( func.function->getRank()==0 || !func.function->hasDerivatives() ) {
    action->error("should have one grid as input to this action");
  }
  // Get the input grid
  func.gridact = ActionWithGrid::getInputActionWithGrid( func.function->getPntrToAction() );
  plumed_assert( func.gridact );
  if( func.gridact->getGridCoordinatesObject().getGridType()!="flat" ) {
    action->error("cannot interpolate on fibonacci sphere");
  }
  std::vector<std::string> argn;
  if( action->keywords.exists("ZERO_OUTSIDE_GRID_RANGE") ) {
    action->parseFlag("ZERO_OUTSIDE_GRID_RANGE",func.set_zero_outside_range);
  } else {
    func.set_zero_outside_range=false;
  }
  if( func.set_zero_outside_range ) {
    action->log.printf("  function is zero outside grid range \n");
  }
  // Get the type of interpolation that we are doing
  std::string itype;
  action->parse("INTERPOLATION_TYPE",itype);
  if( itype=="spline" ) {
    func.interpolation_type=spline;
    func.spline_interpolator=Tools::make_unique<Interpolator>( func.function, func.gridact->getGridCoordinatesObject() );
  } else if( itype=="linear" ) {
    func.interpolation_type=linear;
  } else if( itype=="floor" ) {
    func.interpolation_type=floor;
  } else if( itype=="ceiling" ) {
    func.interpolation_type=ceiling;
  } else {
    action->error("type " + itype + " of interpolation is not defined");
  }
  action->log.printf("  generating off grid points using %s interpolation \n", itype.c_str() );
}

void EvaluateGridFunction::calc( const EvaluateGridFunction& func, bool noderiv, const View<const double> args, function::FunctionOutput& funcout ) {
  const GridCoordinatesObject & gridobject = func.getGridObject();
  if( func.set_zero_outside_range && !gridobject.inbounds( args ) ) {
    funcout.values[0]=0.0;
    return;
  }
  unsigned dimension = gridobject.getDimension();
  plumed_dbg_assert( args.size()==dimension && funcout.values.size()==1 );
  if( func.interpolation_type==spline ) {
    std::vector<double> der( dimension );
    funcout.values[0] =  func.spline_interpolator->splineInterpolation( args, der );
    if( noderiv ) {
      return ;
    }
    for(unsigned j=0; j<dimension; ++j) {
      funcout.derivs[0][j] = der[j];
    }
  } else if( func.interpolation_type==linear ) {
    std::vector<double> xfloor(dimension);
    std::vector<unsigned> indices(dimension), nindices(dimension), ind(dimension);
    gridobject.getIndices( args, indices );
    unsigned nn=gridobject.getIndex(args);
    gridobject.getGridPointCoordinates( nn, nindices, xfloor );
    double y1 = func.function->get(nn);
    funcout.values[0] = y1;
    for(unsigned i=0; i<args.size(); ++i) {
      int x0=1;
      if(nindices[i]==indices[i]) {
        x0=0;
      }
      double ddx=gridobject.getGridSpacing()[i];
      double X = fabs((args[i]-xfloor[i])/ddx-(double)x0);
      for(unsigned j=0; j<args.size(); ++j) {
        ind[j] = indices[j];
      }
      if( gridobject.isPeriodic(i) && (ind[i]+1)==gridobject.getNbin(false)[i] ) {
        ind[i]=0;
      } else {
        ind[i] = ind[i] + 1;
      }
      if( X<epsilon ) {
        continue ;
      }
      funcout.values[0] += ( func.function->get( gridobject.getIndex(ind) ) - y1 )*X;
      if( noderiv ) {
        continue ;
      }
      funcout.derivs[0][i] = ( func.function->get( gridobject.getIndex(ind) ) - y1 ) / ddx;
    }
  } else if( func.interpolation_type==floor ) {
    std::vector<unsigned> indices(dimension);
    gridobject.getIndices( args, indices );
    unsigned nn = gridobject.getIndex(indices);
    plumed_dbg_assert( nn<func.function->getNumberOfValues() );
    funcout.values[0] = func.function->get( nn );
    if( noderiv ) {
      return;
    }
    for(unsigned j=0; j<dimension; ++j) {
      funcout.derivs[0][j] = func.function->getGridDerivative( nn, j );
    }
  } else if( func.interpolation_type==ceiling ) {
    std::vector<unsigned> indices(dimension);
    gridobject.getIndices( args, indices );
    for(unsigned i=0; i<indices.size(); ++i) {
      if( gridobject.isPeriodic(i) && (indices[i]+1)==gridobject.getNbin(false)[i] ) {
        indices[i]=0;
      } else {
        indices[i] = indices[i] + 1;
      }
    }
    unsigned nn = gridobject.getIndex(indices);
    funcout.values[0] = func.function->get( nn );
    if( noderiv ) {
      return;
    }
    for(unsigned j=0; j<dimension; ++j) {
      funcout.derivs[0][j] = func.function->getGridDerivative( nn, j );
    }
  } else {
    plumed_error();
  }
}

}
}
