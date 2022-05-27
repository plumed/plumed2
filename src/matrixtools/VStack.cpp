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
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace matrixtools {

class VStack :
  public ActionWithArguments,
  public ActionWithValue
{
private: 
  std::vector<double> forcesToApply;
  std::vector<std::string> actionsLabelsInChain;
public:
  static void registerKeywords( Keywords& keys );
  explicit VStack(const ActionOptions&);
  unsigned getNumberOfDerivatives() const override;
  unsigned getNumberOfColumns() const override;
  unsigned getNumberOfFinalTasks() override;
  std::vector<unsigned> getMatrixShapeForFinalTasks() override; 
  void calculate() override; 
  void update() override;
  void runFinalJobs() override;
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
  bool performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const ;
  void apply() override;
};

PLUMED_REGISTER_ACTION(VStack,"VSTACK")

void VStack::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); 
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG"); 
}

VStack::VStack(const ActionOptions&ao):
Action(ao),
ActionWithArguments(ao),
ActionWithValue(ao)
{
   std::vector<unsigned> shape( getMatrixShapeForFinalTasks() ); 
   std::vector<Value*> args( getArguments() ); requestArguments( args, true );
   addValue( shape ); bool periodic=false; std::string smin, smax; 
   if( getPntrToArgument(0)->isPeriodic() ) { 
       periodic=true; getPntrToArgument(0)->getDomain( smin, smax );
       setPeriodic( smin, smax );
   } else setNotPeriodic(); 

   for(unsigned i=0; i<getNumberOfArguments();++i) {
      if( periodic && !getPntrToArgument(i)->isPeriodic() ) error("one argument is periodic but " + getPntrToArgument(i)->getName() + " is not periodic");
      if( !periodic && getPntrToArgument(i)->isPeriodic() ) error("one argument is not periodic but " + getPntrToArgument(i)->getName() + " is periodic");
      if( periodic && getPntrToArgument(i)->isPeriodic() ) {
          std::string tmin, tmax; getPntrToArgument(i)->getDomain( tmin, tmax );
          if( tmin!=smin || tmax!=smax ) error("domain of argument " + getPntrToArgument(i)->getName() + " is different from domain for all other arguments");
      }
  }
  forcesToApply.resize( shape[0]*shape[1] ); 
  for(unsigned i=0;i<getNumberOfFinalTasks();++i) addTaskToList(i);

  bool usechain=(shape.size()==2 && distinct_arguments.size()>0);
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      // This ensures we do not try to chain off collect frames 
      if( getPntrToArgument(i)->isTimeSeries() ) { usechain=false; getPntrToComponent(0)->makeHistoryDependent(); break; }
  }
  if( usechain ) { unsigned nd = setupActionInChain(0); forcesToApply.resize(nd); }
}

unsigned VStack::getNumberOfDerivatives() const {
  return forcesToApply.size();
}

unsigned VStack::getNumberOfColumns() const {
  return getPntrToOutput(0)->getShape()[1];
}

unsigned VStack::getNumberOfFinalTasks() {
  std::vector<unsigned> shape( getMatrixShapeForFinalTasks() );
  return shape[0];
}

std::vector<unsigned> VStack::getMatrixShapeForFinalTasks() {
  if( getNumberOfArguments()==0 ) error("no arguments have been selected");
  unsigned nvals = getPntrToArgument(0)->getNumberOfValues();;
  // Check for consistent numbers of values in other actions
  for(unsigned j=0;j<getNumberOfArguments();++j) {
      if( getPntrToArgument(j)->getNumberOfValues()!=nvals ) error("mismatch between number of values in each vector that is to be combined");
  }

  std::vector<unsigned> shape(2);
  if( getPntrToArgument(0)->getRank()==0 ) { shape[0]=1; shape[1]=getNumberOfArguments(); }
  else { shape[0]=nvals; shape[1]=getNumberOfArguments(); }
  return shape; 
}

void VStack::runFinalJobs() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  resizeForFinalTasks();
  runAllTasks();
}

void VStack::update() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  if( getFullNumberOfTasks()>0 ) runAllTasks();
}

void VStack::calculate() {
  // Everything is done elsewhere
  if( actionInChain() || skipCalculate() ) return;
  // Run all the tasks
  runAllTasks();
}

bool VStack::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  // Do not perform the loop here with a loop over other matrix elements
  if( controller!=getLabel() ) return false;

  unsigned iarg = index2 - getFullNumberOfTasks();
  unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
  unsigned istrn = getPntrToArgument( iarg )->getPositionInStream(); 
  if( actionInChain() ) myvals.addValue( ostrn, myvals.get( istrn ) );
  else myvals.addValue( ostrn, getPntrToArgument(iarg)->get( index1 ) );  

  if( doNotCalculateDerivatives() ) return true;

  if( actionInChain() ) { 
      for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
          unsigned kind=myvals.getActiveIndex(istrn,k);
          myvals.addDerivative( ostrn, arg_deriv_starts[iarg] + kind, myvals.getDerivative( istrn, kind ) );
          myvals.updateIndex( ostrn, arg_deriv_starts[iarg] + kind );
      }
  } else if( getPntrToArgument(iarg)->getRank()==0 ) {
      myvals.addDerivative( ostrn, iarg, 1 ); myvals.updateIndex( ostrn, iarg );
  } else {
      unsigned matind = iarg*getPntrToArgument(0)->getNumberOfValues() + index1;
      myvals.addDerivative( ostrn, matind, 1 ); myvals.updateIndex( ostrn, matind );
  }

  return true;
}

void VStack::performTask( const unsigned& current, MultiValue& myvals ) const {
  unsigned base = getFullNumberOfTasks();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    // This does everything in the stream that is done with single matrix elements
    runTask( getLabel(), myvals.getTaskIndex(), current, base + i, myvals );
    // Now clear only elements that are not accumulated over whole row 
    clearMatrixElements( myvals );
  } 
  if( doNotCalculateDerivatives() ) return;
  // Now update the matrix indices
  unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash();
  std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( nmat ) ); unsigned ntot_mat=0;
  if( mat_indices.size()<getNumberOfDerivatives() ) mat_indices.resize( getNumberOfDerivatives() );
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()==0 ) {
          ntot_mat++; mat_indices[ntot_mat] = i;
      } else if( !actionInChain() ) {
          mat_indices[ntot_mat] = i*getPntrToArgument(0)->getNumberOfValues() + current; ntot_mat++;
      } else {
          // Ensure we only store one lot of derivative indices
          bool found=false;
          for(unsigned j=0; j<i; ++j) {
              if( arg_deriv_starts[j]==arg_deriv_starts[i] ) { found=true; break; }
          } 
          if( found ) continue;
          unsigned istrn = getPntrToArgument(i)->getPositionInStream();
          for(unsigned k=0; k<myvals.getNumberActive( istrn ); ++k) mat_indices[ntot_mat + k] = arg_deriv_starts[i] + myvals.getActiveIndex(istrn,k);
          ntot_mat += myvals.getNumberActive( istrn );
      }
  }
  myvals.setNumberOfMatrixIndices( nmat, ntot_mat );
}

void VStack::apply() {
  if( doNotCalculateDerivatives() ) return;

  if( forcesToApply.size()!=getNumberOfDerivatives() ) forcesToApply.resize( getNumberOfDerivatives() );
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, mm );
}

}
}
