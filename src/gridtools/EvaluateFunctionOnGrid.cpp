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
#include "core/PlumedMain.h"
#include "ActionWithInputGrid.h"

namespace PLMD {
namespace gridtools {

class EvaluateFunctionOnGrid : public ActionWithInputGrid {
private:
  unsigned nderivatives;
  std::vector<double> forcesToApply;
public:
  static void registerKeywords( Keywords& keys );
  explicit EvaluateFunctionOnGrid(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const ;
  void prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList );
  void finishOutputSetup(){}
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
  void apply();
};

PLUMED_REGISTER_ACTION(EvaluateFunctionOnGrid,"EVALUATE_FUNCTION_FROM_GRID")

void EvaluateFunctionOnGrid::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.addFlag("ZERO_OUTSIDE_GRID_RANGE",false,"if we are asked to evaluate the function for a number that is outside the range of the grid set it to zero");
}

EvaluateFunctionOnGrid::EvaluateFunctionOnGrid(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao)
{ 
  Value* gval=getPntrToArgument(0); ActionWithValue* act=gval->getPntrToAction();
  std::vector<unsigned> ind( gval->getRank() ), nbin( gval->getRank() ); std::string gtype;
  std::vector<double> spacing( gval->getRank() ), xx( gval->getRank() ); std::vector<bool> pbc( gval->getRank() );
  std::vector<std::string> argn( gval->getRank() ), min( gval->getRank() ), max( gval->getRank() );
  act->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, false );
  if( gtype=="fibonacci" ) error("cannot interpolate on fibonacci sphere");
  // Arg ends must be setup once more
  arg_ends.resize(0); arg_ends.push_back(1); 
  for(unsigned i=0;i<gval->getRank();++i) arg_ends.push_back( 2+i );
  // Retreive values with these arguments
  std::vector<Value*> argv; ActionWithArguments::interpretArgumentList( argn, plumed.getActionSet(), this, argv ); 
  // Now check that arguments make sense
  for(unsigned i=1;i<argv.size();++i) {
      if( argv[0]->getRank()!=argv[i]->getRank() ) error("mismatched ranks for arguments");
      for(unsigned j=0;j<argv[0]->getRank();++j) {
          if( argv[0]->getShape()[j]!=argv[i]->getShape()[j] ) error("mismatched shapes for arguments");
      } 
  }
  // Now re-request the arguments 
  std::vector<Value*> req_arg; req_arg.push_back( getPntrToArgument(0) );  
  log.printf("  arguments for grid are");
  for(unsigned i=0;i<argv.size();++i) {
      log.printf(" %s", argv[i]->getName().c_str() );  
      req_arg.push_back( argv[i] );
  }
  log.printf("\n"); requestArguments( req_arg, true, 1 );
  parseFlag("ZERO_OUTSIDE_GRID_RANGE",set_zero_outside_range);
  if( set_zero_outside_range ) log.printf("  function is zero outside grid range \n");

  // Make a task list
  createTasksFromArguments();
  // Get the number of derivatives
  nderivatives = getNumberOfArguments() - 1;
  if( distinct_arguments.size()>0 ) nderivatives = setupActionInChain(1);
  forcesToApply.resize( nderivatives );
  // Create value for this function
  std::vector<unsigned> shape(0);
  if( argv[0]->getRank()>0 ) {
      shape.resize( argv[0]->getRank() );
      for(unsigned i=0;i<shape.size();++i) shape[i]=argv[0]->getShape()[i];
      addValue( shape );
  } else addValueWithDerivatives( shape ); 
  setNotPeriodic();
}

unsigned EvaluateFunctionOnGrid::getNumberOfDerivatives() const {
  return nderivatives;
}

void EvaluateFunctionOnGrid::prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList ) {
  if( firststep ) setupGridObject();
}

void EvaluateFunctionOnGrid::performTask( const unsigned& current, MultiValue& myvals ) const {
  unsigned nargs = getNumberOfArguments() - 1, ostrn = getPntrToOutput(0)->getPositionInStream();
  std::vector<double> args( nargs ), der( nargs ); retrieveArguments( myvals, args, 1 );
  // Evaluate function
  double func = getFunctionValueAndDerivatives( args, der ); myvals.addValue( ostrn, func );
  // And derivatives
  if( doNotCalculateDerivatives() ) return;
  // Set the final derivatives
  if( actionInChain() ) {
    for(unsigned j=0; j<nargs; ++j) {
      unsigned istrn = getArgumentPositionInStream(j+1, myvals);
      for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
        unsigned kind=myvals.getActiveIndex(istrn,k);
        myvals.addDerivative( ostrn, arg_deriv_starts[j] + kind, der[j]*myvals.getDerivative( istrn, kind ) );
        myvals.updateIndex( ostrn, arg_deriv_starts[j] + kind );
      }
    }
    return;
  } 
  for(unsigned i=0;i<nargs;++i) { 
      if( getPntrToArgument(i+1)->getRank()>0 ) error("non-chained function has not been implemented here");
      myvals.addDerivative( ostrn, i, der[i] ); myvals.updateIndex( ostrn, i ); 
  }
}

void EvaluateFunctionOnGrid::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned fstart=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 1, forcesToApply, fstart );
}

}
}
