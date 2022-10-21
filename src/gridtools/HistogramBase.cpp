/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "HistogramBase.h"
#include "core/PlumedMain.h"
#include "core/ActionShortcut.h"

namespace PLMD {
namespace gridtools {

void HistogramBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("optional","HEIGHTS","this keyword takes the label of an action that calculates a vector of values.  The elements of this vector "
           "are used as weights for the Gaussians.");
  keys.addFlag("UNORMALIZED",false,"calculate the unormalized distribution of colvars");
}

HistogramBase::HistogramBase(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  hasheight(false),
  fixed_width(false),
  grid_dimension(getNumberOfArguments()),
  numberOfKernels(1)
{
  // Check all the values have the right size
  setNumberOfKernels();
  // Get the heights if need be
  std::vector<std::string> weight_str; parseVector("HEIGHTS",weight_str);
  if( weight_str.size()==1 ) {
    std::vector<Value*> weight_args; ActionWithArguments::interpretArgumentList( weight_str, plumed.getActionSet(), this, weight_args );
    hasheight=true; std::vector<Value*> args( getArguments() ); args.push_back( weight_args[0] );
    log.printf("  quantities used for weights are : %s \n", weight_str[0].c_str() );

    if( weight_args[0]->getNumberOfValues()>1 && numberOfKernels!=weight_args[0]->getNumberOfValues() ) error("mismatch between numbers of values in input arguments and HEIGHTS");
    requestArguments( args, true );
  }
  // Make sure we are storing all the values
  done_over_stream = false;
  for(unsigned i=0;i<getNumberOfArguments();++i) getPntrToArgument(i)->buildDataStore( getLabel() );

  parseFlag("UNORMALIZED",unorm);
  if( unorm ) log.printf("  calculating unormalized distribution \n");
  else log.printf("  calculating normalized distribution \n");
  resizeForcesToApply();
}

void HistogramBase::resizeForcesToApply() {
  // Resize the forces vector
  unsigned nvals_t=0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) nvals_t += getPntrToArgument(i)->getNumberOfValues();
  forcesToApply.resize( nvals_t );
}

void HistogramBase::setNumberOfKernels() {
  numberOfKernels=getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=1; i<grid_dimension; ++i) {
    if( numberOfKernels!=getPntrToArgument(i)->getNumberOfValues() ) error("mismatch between numbers of values in input arguments");
  }
  if( numberOfKernels>1 ) {
      task_list.clear(); bool symmetric=false;
      if( getPntrToArgument(0)->getRank()==2 ) {
          symmetric=getPntrToArgument(0)->isSymmetric();
          for(unsigned i=0; i<grid_dimension; ++i) {
              if( !getPntrToArgument(i)->isSymmetric() ) symmetric=false;
          }
      }
      if( symmetric ) { 
          std::vector<unsigned> shape( getPntrToArgument(0)->getShape() );
          for(unsigned j=1;j<shape[0];++j) {
              for(unsigned k=0;k<j;++k) task_list.insert(AtomNumber::index(j*shape[0]+k));
          } 
      } else {
          for(unsigned i=0;i<numberOfKernels;++i) task_list.insert(AtomNumber::index(i));
      }
  }
}

void HistogramBase::addValueWithDerivatives( const std::vector<unsigned>& shape ) {
  ActionWithValue::addValueWithDerivatives( shape ); setNotPeriodic();
  getPntrToOutput(0)->alwaysStoreValues(); 
  getPntrToOutput(0)->setDerivativeIsZeroWhenValueIsZero();
}

unsigned HistogramBase::getNumberOfDerivatives() const {
  return grid_dimension;
}

void HistogramBase::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  gridobject.getGridPointCoordinates( ind, indices, coords );
}

void HistogramBase::getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const {
  if( setlength ) gridobject.putCoordinateAtValue( ind, getPntrToOutput(0)->get(ind), coords );
  else  gridobject.putCoordinateAtValue( ind, 1.0, coords );
}

void HistogramBase::calculate() {
  // Everything is done elsewhere
  if( actionInChain() ) return;
  // This is done if we are calculating a function of multiple cvs
  runAllTasks();
}

void HistogramBase::update() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  // This is done if we are doing a histogram from a time series
  // Make sure tasks are set
  runAllTasks();
}

void HistogramBase::runFinalJobs() {
  if( skipUpdate() ) return;
  plumed_assert( !actionInChain() );
  // Need to create tasks here
  if( getPntrToArgument(0)->getRank()>0 ) { 
      setNumberOfKernels(); 
      // Make sure tasks are set
      setupCurrentTaskList(); runAllTasks();
  }
} 

void HistogramBase::setupCurrentTaskList() {
  if( !fixed_width ) setupNeighborsVector();
  if( numberOfKernels>1 ) {
    unsigned hind = getNumberOfDerivatives();
    getPntrToOutput(0)->setNumberOfTasks( numberOfKernels );
    if( hasheight && getPntrToArgument(grid_dimension)->getRank()>0 ) {
        for(const auto & t : task_list ) {
            if( fabs(getPntrToArgument(grid_dimension)->get(t.index()))>epsilon ) getPntrToOutput(0)->addTaskToCurrentList(t);
        }
    } else { getPntrToOutput(0)->addTasksToCurrentList( task_list ); }
    norm = static_cast<double>( numberOfKernels );
  } else {
    std::vector<double> args( getNumberOfDerivatives() );
    double height=1.0; if( hasheight ) height = getPntrToArgument(grid_dimension)->get();
    for(unsigned i=0; i<args.size(); ++i) args[i]=getPntrToArgument(i)->get();
    buildSingleKernel( height, args );
  }
}

void HistogramBase::performTask( const unsigned& current, MultiValue& myvals ) const {
  if( numberOfKernels==1 ) {
    std::vector<double> args( getNumberOfDerivatives() ), der( getNumberOfDerivatives() );
    unsigned valout = getPntrToOutput(0)->getPositionInStream();
    gridobject.getGridPointCoordinates( current, args ); double vv = calculateValueOfSingleKernel( args, der );
    myvals.setValue( valout, vv );
    for(unsigned i=0; i<der.size(); ++i) { myvals.addDerivative( valout, i, der[i] ); myvals.updateIndex( valout, i ); }
  }
}

void HistogramBase::retrieveArgumentsAndHeight( const MultiValue& myvals, std::vector<double>& args, double& height ) const {
  height=1.0; for(unsigned i=0; i<args.size(); ++i) args[i]=getPntrToArgument(i)->get( myvals.getTaskIndex() );
  if( hasheight && getPntrToArgument(grid_dimension)->getRank()==0 ) height = getPntrToArgument( grid_dimension )->get();
  else if( hasheight ) height = height = getPntrToArgument( grid_dimension )->get( myvals.getTaskIndex() );
  if( !unorm ) height = height / norm;
}

void HistogramBase::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                       const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( valindex==0 );
  if( numberOfKernels==1 ) {
    unsigned istart = bufstart + (1+getNumberOfDerivatives())*code;
    unsigned valout = getPntrToOutput(0)->getPositionInStream(); buffer[istart] += myvals.get( valout );
    for(unsigned i=0; i<getNumberOfDerivatives(); ++i) buffer[istart+1+i] += myvals.getDerivative( valout, i );
    return;
  }

  // Add the kernel to the grid
  std::vector<double> args( getNumberOfDerivatives() ); double height;
  retrieveArgumentsAndHeight( myvals, args, height );
  if( fabs(height)>epsilon ) addKernelToGrid( height, args, bufstart, buffer );
}

void HistogramBase::apply() {
  // Everything is done elsewhere
  if( doNotCalculateDerivatives() ) return;
  // And add forces
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned ss=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, ss );
}

void HistogramBase::gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  if( numberOfKernels==1 ) {
    if( getPntrToOutput(0)->forcesWereAdded() ) {
      unsigned valout = getPntrToOutput(0)->getPositionInStream(); double fforce = getPntrToOutput(0)->getForce( itask );
      for(unsigned i=0; i<getNumberOfDerivatives(); ++i) forces[i] += fforce*myvals.getDerivative( valout, i );
    }
    return;
  }
  std::vector<double> args( getNumberOfDerivatives() ); double height;
  retrieveArgumentsAndHeight( myvals, args, height );
  if( fabs(height)>epsilon ) {
    if( getPntrToOutput(0)->forcesWereAdded() ) addKernelForces( hasheight, itask, args, itask, height, forces );
  }
}

}
}
