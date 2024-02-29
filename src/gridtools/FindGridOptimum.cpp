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

//+PLUMEDOC GRIDCALC FIND_GRID_MINIMUM
/*
Find the point with the lowest value of the function on the grid

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC GRIDCALC FIND_GRID_MAXIMUM
/*
Find the point with the highest value of the function on the grid

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class FindGridOptimum : public ActionWithGrid {
private:
  bool domin;
  double cgtol;
  std::unique_ptr<Interpolator> function;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindGridOptimum(const ActionOptions&ao);
  void setupOnFirstStep( const bool incalc ) override { plumed_error(); }
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  std::vector<std::string> getGridCoordinateNames() const override ;
  const GridCoordinatesObject& getGridCoordinatesObject() const override ;
  double calculateValueAndDerivatives( const std::vector<double>& pp, std::vector<double>& der );
  void performTask( const unsigned& current, MultiValue& myvals ) const override { plumed_error(); }
};

PLUMED_REGISTER_ACTION(FindGridOptimum,"FIND_GRID_MAXIMUM")
PLUMED_REGISTER_ACTION(FindGridOptimum,"FIND_GRID_MINIMUM")

void FindGridOptimum::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords( keys ); keys.use("ARG");
  keys.addFlag("NOINTERPOL",false,"do not interpolate the function when finding the optimum");
  keys.add("compulsory","CGTOL","1E-4","the tolerance for the conjugate gradient optimization");
  keys.addOutputComponent("optval","default","the value of the function at the optimum");
  keys.addOutputComponent("_opt","default","the values of the arguments of the function at the optimum can be referenced elsewhere in the input file "
                          "by using the names of the arguments followed by the string _opt");
}

FindGridOptimum::FindGridOptimum(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  cgtol(0)
{
  if( getName()=="FIND_GRID_MAXIMUM" ) domin=false;
  else if( getName()=="FIND_GRID_MINIMUM" ) domin=true;
  else plumed_error();
  // Create value for this function
  std::vector<std::string> argn( getGridCoordinateNames() );
  std::vector<unsigned> shape(0); for(unsigned i=0; i<argn.size(); ++i) addComponent( argn[i] + "_opt", shape );
  addComponent( "optval", shape ); componentIsNotPeriodic( "optval" );
  bool nointerpol=false; parseFlag("NOINTERPOL",nointerpol);
  if( !nointerpol ) { parse("CGTOL",cgtol); function=Tools::make_unique<Interpolator>( getPntrToArgument(0), getGridCoordinatesObject() ); }
}

unsigned FindGridOptimum::getNumberOfDerivatives() {
  return 0;
}

const GridCoordinatesObject& FindGridOptimum::getGridCoordinatesObject() const {
  ActionWithGrid* ag=ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag ); return ag->getGridCoordinatesObject();
}

std::vector<std::string> FindGridOptimum::getGridCoordinateNames() const {
  ActionWithGrid* ag=ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag ); return ag->getGridCoordinateNames();
}

double FindGridOptimum::calculateValueAndDerivatives( const std::vector<double>& pp, std::vector<double>& der ) {
  double val = function->splineInterpolation( pp, der );
  // We normalise the derivatives here and set them so that the linesearch is done over the cell that we know
  // in the grid that we know the minimum is inside
  std::vector<double> spacing( getGridCoordinatesObject().getGridSpacing() );
  double norm = 0; for(unsigned i=0; i<der.size(); ++i) norm += der[i]*der[i];
  norm = sqrt(norm); for(unsigned i=0; i<der.size(); ++i) der[i] = spacing[i]*der[i] / norm;
  if( domin ) return val;
  // If we are looking for a maximum
  for(unsigned i=0; i<der.size(); ++i) der[i] = -der[i];
  return -val;
}

void FindGridOptimum::calculate() {
  const GridCoordinatesObject& ingrid = getGridCoordinatesObject(); Value* gval = getPntrToArgument(0);
  std::vector<double> optargs( gval->getRank() ); std::vector<unsigned> gridind( gval->getRank() );
  double optval=gval->get( 0 ); ingrid.getGridPointCoordinates( 0, gridind, optargs );
  unsigned nval = gval->getNumberOfValues(); bool constant=true;
  for(unsigned i=0; i<nval; ++i) {
    double tval = gval->get( i );
    if( domin && (tval<optval || std::isnan(optval)) ) { constant=false; optval=tval; ingrid.getGridPointCoordinates( i, gridind, optargs ); }
    if( !domin && (tval>optval || std::isnan(optval)) ) { constant=false; optval=tval; ingrid.getGridPointCoordinates( i, gridind, optargs ); }
  }
  // This basically ensures we deal with cases where all points on the grid are infinity as isinf doesn't work on intel compiler
  if( constant ) {
    if( domin && gval->get(0)>=gval->get(1) ) return;
    else if( gval->get(0)<=gval->get(1) ) return;
  }
  if( std::isinf(optval) ) { return; }

  if( std::isnan(optval) ) error("all values on grid are nans");
  // And do conjugate gradient optimisation (because we can!!)
  if( cgtol>0 ) {
    ConjugateGradient<FindGridOptimum> myminimiser( this );
    myminimiser.minimise( cgtol, optargs, &FindGridOptimum::calculateValueAndDerivatives );
  }
  // And set the final value
  for(unsigned j=0; j<optargs.size(); ++j) getPntrToComponent(j)->set( optargs[j] );
  std::vector<double> optder( gval->getRank() );
  getPntrToComponent(optargs.size())->set( function->splineInterpolation( optargs, optder ) );
}

}
}
