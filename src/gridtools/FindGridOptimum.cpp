/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "tools/ConjugateGradient.h"
#include "ActionWithGrid.h"
#include "EvaluateGridFunction.h"
#include "Interpolator.h"

//+PLUMEDOC GRIDANALYSIS FIND_GRID_MINIMUM
/*
Find the point with the lowest value of the function on the grid

This command takes a function on a grid as input and returns multiple scalars.  One of the returned scalars (the one with component `optval`) is the value of the input
function at its lowest point. The other scalars returned are the coordinates for which the function takes this particular value.

As the input grid has the values of the function evaluated over a grid of points we find this minimum value by finding the point on the grid where the function has its lowest value using
interpolation and conjugate gradients.  The coordinates of the grid point where the function's value is lowest will be used as an initial coordinate
for the optimisation.  We then interpolate the function to find the optimum.  You can use use the CGTOL keyword to control the conjugate gradient algorithm that is used to find the location of the
minimum in the interpolated function.

To print out the location on the grid where the function is minimised and the value of the function at that point you can use an input like the one shown below:

```plumed
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
ff: CONVERT_TO_FES ARG=hA1 TEMP=300
min: FIND_GRID_MINIMUM ARG=ff
PRINT ARG=min.x_opt,min.optval STRIDE=0
```

Notice that we set STRIDE=0 in the PRINT command here so the position of the minimum is only output once at the end of the simulation.

A more common usage of this action is illustrated in the following input, which demonstrates how you can use FIND_GRID_MINIMUM to shift a free energy surface so that the minimum in the output
fes has a value of zero.

```plumed
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
ff: CONVERT_TO_FES ARG=hA1 TEMP=300
min: FIND_GRID_MINIMUM ARG=ff NOINTERPOL
sff: CUSTOM ARG=ff,min.optval FUNC=x+y PERIODIC=NO
DUMPGRID ARG=sff FILE=fes.dat
```

Notice that the `NOINTERPOL` flag has been used here so the interpolation and conjugate gradient algorithm is not used. The output is simply the grid point at which the function takes its lowest value.

*/
//+ENDPLUMEDOC

//+PLUMEDOC GRIDANALYSIS FIND_GRID_MAXIMUM
/*
Find the point with the highest value of the function on the grid

This command takes a function on a grid as input and returns multiple scalars.  One of the returned scalars (the one with component `optval`) is the value of the input
function at its highest point. The other scalars returned are the coordinates for which the function takes this particular value.

As the input grid has the values of the function evaluated over a grid of points we find this maximum value by finding the point on the grid where the function has its lowest value using
interpolation and conjugate gradients.  The coordinates of the grid point where the function's value is highest will be used as an initial coordinate
for the optimisation.  We then interpolate the function to find the optimum.  You can use use the CGTOL keyword to control the conjugate gradient algorithm that is used to find the location of the
maximum in the interpolated function.

To print out the location on the grid where the function is maximised and the value of the function at that point you can use an input like the one shown below:

```plumed
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
max: FIND_GRID_MAXIMUM ARG=hA1
PRINT ARG=max.x_opt,max.optval STRIDE=0
```

Notice that we set STRIDE=0 in the PRINT command here so the position of the maximum is only output once at the end of the simulation.

If for any reason you do not want to use interpolation and conjugate gradient to find the optimum you can use the `NOINTERPOL` flag as shown below.

```plumed
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
max: FIND_GRID_MAXIMUM ARG=hA1 NOINTERPOL
PRINT ARG=max.x_opt,max.optval STRIDE=0
```

The FIND_GRID_MAXIMUM action then outputs the coordinates of the grid point where the function is largest as well as the value of the function on that grid point.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class FindGridOptimum :
  public ActionWithValue,
  public ActionWithArguments {
private:
  bool domin;
  double cgtol;
  std::vector<double> spacing;
  std::unique_ptr<Interpolator> function;
  double calculateValueAndDerivatives( const std::vector<double>& pp, std::vector<double>& der );
public:
  static void registerKeywords( Keywords& keys );
  explicit FindGridOptimum(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(FindGridOptimum,"FIND_GRID_MAXIMUM")
PLUMED_REGISTER_ACTION(FindGridOptimum,"FIND_GRID_MINIMUM")

void FindGridOptimum::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addInputKeyword("compulsory","ARG","grid","the label for the function on the grid that you would like to find the optimum in");
  keys.addFlag("NOINTERPOL",false,"do not interpolate the function when finding the optimum");
  keys.add("compulsory","CGTOL","1E-4","the tolerance for the conjugate gradient optimization");
  keys.addOutputComponent("optval","default","scalar","the value of the function at the optimum");
  keys.addOutputComponent("_opt","default","scalar","the values of the arguments of the function at the optimum can be referenced elsewhere in the input file "
                          "by using the names of the arguments followed by the string _opt");
}

FindGridOptimum::FindGridOptimum(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  cgtol(0) {
  if( getName()=="FIND_GRID_MAXIMUM" ) {
    domin=false;
  } else if( getName()=="FIND_GRID_MINIMUM" ) {
    domin=true;
  } else {
    plumed_error();
  }
  // Create value for this function
  ActionWithGrid* ag=ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  if( !ag ) {
    error("input action should be a grid");
  }
  std::vector<std::string> argn( ag->getGridCoordinateNames() );
  std::vector<std::size_t> shape(0);
  for(unsigned i=0; i<argn.size(); ++i) {
    addComponent( argn[i] + "_opt", shape );
    componentIsNotPeriodic( argn[i] + "_opt" );
  }
  addComponent( "optval", shape );
  componentIsNotPeriodic( "optval" );
  bool nointerpol=false;
  parseFlag("NOINTERPOL",nointerpol);
  if( !nointerpol ) {
    parse("CGTOL",cgtol);
    function=Tools::make_unique<Interpolator>( getPntrToArgument(0), ag->getGridCoordinatesObject() );
  }
}

unsigned FindGridOptimum::getNumberOfDerivatives() {
  return 0;
}

double FindGridOptimum::calculateValueAndDerivatives( const std::vector<double>& pp, std::vector<double>& der ) {
  double val = function->splineInterpolation( pp, der );
  // We normalise the derivatives here and set them so that the linesearch is done over the cell that we know
  // in the grid that we know the minimum is inside
  double norm = 0;
  for(unsigned i=0; i<der.size(); ++i) {
    norm += der[i]*der[i];
  }
  norm = sqrt(norm);
  for(unsigned i=0; i<der.size(); ++i) {
    der[i] = spacing[i]*der[i] / norm;
  }
  if( domin ) {
    return val;
  }
  // If we are looking for a maximum
  for(unsigned i=0; i<der.size(); ++i) {
    der[i] = -der[i];
  }
  return -val;
}

void FindGridOptimum::calculate() {
  ActionWithGrid* ag=ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag );
  spacing = ag->getGridCoordinatesObject().getGridSpacing();
  const GridCoordinatesObject& ingrid = ag->getGridCoordinatesObject();
  Value* gval = getPntrToArgument(0);
  std::vector<double> optargs( gval->getRank() );
  std::vector<unsigned> gridind( gval->getRank() );
  double optval=gval->get( 0 );
  ingrid.getGridPointCoordinates( 0, gridind, optargs );
  unsigned nval = gval->getNumberOfValues();
  bool constant=true;
  for(unsigned i=0; i<nval; ++i) {
    double tval = gval->get( i );
    if( domin && (tval<optval || std::isnan(optval)) ) {
      constant=false;
      optval=tval;
      ingrid.getGridPointCoordinates( i, gridind, optargs );
    }
    if( !domin && (tval>optval || std::isnan(optval)) ) {
      constant=false;
      optval=tval;
      ingrid.getGridPointCoordinates( i, gridind, optargs );
    }
  }
  // This basically ensures we deal with cases where all points on the grid are infinity as isinf doesn't work on intel compiler
  if( constant ) {
    if( domin && gval->get(0)>=gval->get(1) ) {
      return;
    } else if( gval->get(0)<=gval->get(1) ) {
      return;
    }
  }
  if( std::isinf(optval) ) {
    return;
  }

  if( std::isnan(optval) ) {
    error("all values on grid are nans");
  }
  // And do conjugate gradient optimisation (because we can!!)
  if( cgtol>0 ) {
    ConjugateGradient<FindGridOptimum> myminimiser( this );
    myminimiser.minimise( cgtol, optargs, &FindGridOptimum::calculateValueAndDerivatives );
  }
  // And set the final value
  for(unsigned j=0; j<optargs.size(); ++j) {
    getPntrToComponent(j)->set( optargs[j] );
  }
  std::vector<double> optder( gval->getRank() );
  getPntrToComponent(optargs.size())->set( function->splineInterpolation( optargs, optder ) );
}

}
}
