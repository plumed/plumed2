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
#include "ActionWithInputGrid.h"

namespace PLMD {
namespace gridtools {

class FindGridOptimum : public ActionWithInputGrid {
private: 
  bool domin;
  double cgtol;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindGridOptimum(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const ;
  double calculateValueAndDerivatives( const std::vector<double>& pp, std::vector<double>& der );
  void finishOutputSetup() {}
  void runTheCalculation();
  void apply() {}
};

PLUMED_REGISTER_ACTION(FindGridOptimum,"FIND_GRID_MAXIMUM")
PLUMED_REGISTER_ACTION(FindGridOptimum,"FIND_GRID_MINIMUM")

void FindGridOptimum::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.addFlag("NOINTERPOL",false,"do not interpolate the function when finding the optimum");
  keys.add("compulsory","CGTOL","1E-4","the tolerance for the conjugate gradient optimization");
  keys.addOutputComponent("optval","default","the value of the function at the optimum");
  keys.addOutputComponent("_opt","default","the values of the arguments of the function at the optimum can be referenced elsewhere in the input file "
                                           "by using the names of the arguments followed by the string _opt");
}

FindGridOptimum::FindGridOptimum(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao),
cgtol(0)
{ 
  if( getName()=="FIND_GRID_MAXIMUM" ) domin=false;
  else if( getName()=="FIND_GRID_MINIMUM" ) domin=true;
  else plumed_error(); 
  Value* gval=getPntrToArgument(0); ActionWithValue* act=gval->getPntrToAction();
  std::vector<unsigned> ind( gval->getRank() ), nbin( gval->getRank() ); std::string gtype;
  std::vector<double> spacing( gval->getRank() ), xx( gval->getRank() ); std::vector<bool> pbc( gval->getRank() );
  std::vector<std::string> argn( gval->getRank() ), min( gval->getRank() ), max( gval->getRank() );
  act->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, false );
  if( gtype=="fibonacci" ) error("cannot find minimum on fibonacci sphere");
  // Arg ends must be setup once more
  arg_ends.resize(0); arg_ends.push_back(1); 
  // Create value for this function
  std::vector<unsigned> shape(0); 
  for(unsigned i=0;i<argn.size();++i) {
      addComponent( argn[i] + "_opt", shape ); 
      if( pbc[i] ) componentIsPeriodic( argn[i] + "_opt", min[i], max[i] );
      else componentIsNotPeriodic( argn[i] + "_opt" );
  }
  addComponent( "optval", shape ); componentIsNotPeriodic( "optval" );
  bool nointerpol=false; parseFlag("NOINTERPOL",nointerpol); if( !nointerpol ) parse("CGTOL",cgtol);
}

unsigned FindGridOptimum::getNumberOfDerivatives() const {
  return 0;
}

double FindGridOptimum::calculateValueAndDerivatives( const std::vector<double>& pp, std::vector<double>& der ) {
  double val = getFunctionValueAndDerivatives( pp, der );
  // We normalise the derivatives here and set them so that the linesearch is done over the cell that we know 
  // in the grid that we know the minimum is inside
  double norm = 0; for(unsigned i=0;i<der.size();++i) norm += der[i]*der[i];
  norm = sqrt(norm); for(unsigned i=0;i<der.size();++i) der[i] = getGridObject().getGridSpacing()[i]*der[i] / norm;
  if( domin ) return val;
  // If we are looking for a maximum
  for(unsigned i=0;i<der.size();++i) der[i] = -der[i];
  return -val;
}

void FindGridOptimum::runTheCalculation() {
  Value* gval = getPntrToArgument(0); ActionWithValue* gact=gval->getPntrToAction(); 
  std::vector<double> optargs( gval->getRank() ); std::vector<unsigned> gridind( gval->getRank() ); 
  double optval=getFunctionValue( 0 ); gact->getGridPointIndicesAndCoordinates( 0, gridind, optargs );
  unsigned nval = gval->getNumberOfValues();
  for(unsigned i=0;i<nval;++i) {
      double tval = getFunctionValue( i ); 
      if( domin && (tval<optval || std::isnan(optval)) ) { optval=tval; gact->getGridPointIndicesAndCoordinates( i, gridind, optargs ); }
      if( !domin && (tval>optval || std::isnan(optval)) ) { optval=tval; gact->getGridPointIndicesAndCoordinates( i, gridind, optargs ); }
  }
  if( std::isnan(optval) ) error("all values on grid are nans");
  // And do conjugate gradient optimisation (because we can!!)
  if( cgtol>0 ) {
      ConjugateGradient<FindGridOptimum> myminimiser( this );
      myminimiser.minimise( cgtol, optargs, &FindGridOptimum::calculateValueAndDerivatives );
  }
  // And set the final value
  for(unsigned j=0;j<optargs.size();++j) getPntrToOutput(j)->set( optargs[j] );
  std::vector<double> optder( gval->getRank() );
  getPntrToOutput(optargs.size())->set( getFunctionValueAndDerivatives( optargs, optder ) );
}

}
}
