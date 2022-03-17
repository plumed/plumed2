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

void HistogramBase::createKDEObject( const std::string& lab, const std::string& command, const std::string& height, const std::string& height_str, ActionShortcut* action ) {
  std::string inp; bool uflag; action->parseFlag("UNORMALIZED",uflag);
  // Deal with the weights if we are doing averages on a grid
  if( height.length()>0 && !uflag ) {
    inp = lab + "_unorm: " + command + "_CALC " + action->convertInputLineToString(); 
    action->readInputLine( lab + "_hsum: SUM ARG=" + height + " PERIODIC=NO");
    inp = inp + " UNORMALIZED";
  } else if( !uflag ) {
     inp = lab + ": " + command + "_CALC " + action->convertInputLineToString();
  } else {
     inp = lab + ": " + command + "_CALC UNORMALIZED " + action->convertInputLineToString();
  }
  action->readInputLine( inp + height_str );
  if( height.length()>0 && !uflag ) {
    action->readInputLine(  lab + ": MATHEVAL ARG1=" + lab + "_unorm ARG2=" + lab + "_hsum FUNC=x/y PERIODIC=NO"); 
  }
}

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
  heights_index(1),
  numberOfKernels(1),
  one_kernel_at_a_time(false)
{
  // Check all the values have the right size
  if( arg_ends.size()>0 ) {
    setNumberOfKernels();
  } else {
    arg_ends.push_back(0); for(unsigned i=0; i<getNumberOfArguments(); ++i) arg_ends.push_back(i+1);
  }
  // Get the heights if need be
  std::vector<std::string> weight_str; parseVector("HEIGHTS",weight_str);
  if( weight_str.size()>0 ) {
    std::vector<Value*> weight_args; ActionWithArguments::interpretArgumentList( weight_str, plumed.getActionSet(), this, weight_args );
    heights_index=2; std::vector<Value*> args( getArguments() ); unsigned tvals=0;
    log.printf("  quantities used for weights are : %s ", weight_str[0].c_str() );
    for(unsigned i=1; i<weight_args.size(); ++i) log.printf(", %s", weight_str[i].c_str() );
    log.printf("\n");

    for(unsigned i=0; i<weight_args.size(); ++i) {
      tvals += weight_args[i]->getNumberOfValues();
      args.push_back( weight_args[i] );
    }
    if( tvals>1 && numberOfKernels!=tvals ) error("mismatch between numbers of values in input arguments and HEIGHTS");
    arg_ends.push_back( args.size() ); requestArguments( args, true );
  }

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

void HistogramBase::createTaskList() {
  bool hasrank = getPntrToArgument(0)->getRank()>0; done_over_stream = false;
  if( hasrank ) {
    buildTasksFromBasedOnRankOfInputData();
  } else {
    one_kernel_at_a_time=true; for(unsigned i=0; i<arg_ends.size(); ++i) { if( arg_ends[i]!=i ) { one_kernel_at_a_time=false; break; } }
    if( !one_kernel_at_a_time ) for(unsigned i=0; i<numberOfKernels; ++i) addTaskToList(i);
  }
  for(unsigned i=0; i<getNumberOfArguments(); ++i) { getPntrToArgument(i)->buildDataStore( getLabel() ); }
}

void HistogramBase::setNumberOfKernels() {
  numberOfKernels=0; for(unsigned i=arg_ends[0]; i<arg_ends[1]; ++i) numberOfKernels += getPntrToArgument(i)->getNumberOfValues();
  for(unsigned i=1; i<arg_ends.size()-1; ++i) {
    unsigned tvals=0; for(unsigned j=arg_ends[i]; j<arg_ends[i+1]; ++j) tvals += getPntrToArgument(j)->getNumberOfValues();
    if( numberOfKernels!=tvals ) error("mismatch between numbers of values in input arguments");
  }
}

void HistogramBase::buildTasksFromBasedOnRankOfInputData() {
  plumed_dbg_assert( getPntrToArgument(0)->getRank()>0 );
  // Now build the data task list based on the rank of the input data
  if( getPntrToArgument(0)->getRank()==2 ) {
     std::vector<unsigned> shape( getPntrToArgument(0)->getShape() );
     bool symmetric=getPntrToArgument(arg_ends[getNumberOfDerivatives()])->isSymmetric(); 
     for(unsigned i=0; i<getNumberOfDerivatives(); ++i) {
         if( !getPntrToArgument(i)->isSymmetric() ) symmetric=false;
     }
     if( symmetric ) {
         for(unsigned j=0;j<shape[0];++j) {
             for(unsigned k=0;k<=j;++k) addTaskToList( j*shape[0] + k );
         }
     } else {
         for(unsigned i=0; i<numberOfKernels; ++i) addTaskToList(i);
     }
  } else if( getPntrToArgument(0)->getRank()==1 ) {
     for(unsigned i=0; i<numberOfKernels; ++i) addTaskToList(i);
  } else {
     error("do not know how to build histograms for objects with this rank");
  }
}

void HistogramBase::addValueWithDerivatives( const std::vector<unsigned>& shape ) {
  ActionWithValue::addValueWithDerivatives( shape ); setNotPeriodic();
  getPntrToOutput(0)->alwaysStoreValues();
  if( one_kernel_at_a_time ) {
    for(unsigned i=0; i<gridobject.getNumberOfPoints(); ++i) addTaskToList(i);
  }
}

unsigned HistogramBase::getNumberOfDerivatives() const {
  return arg_ends.size()-heights_index;
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
  runAllTasks();
}

void HistogramBase::runFinalJobs() {
  if( skipUpdate() ) return;
  plumed_assert( !actionInChain() );
  // Need to create tasks here
  if( getFullNumberOfTasks()==0 ) { setNumberOfKernels(); buildTasksFromBasedOnRankOfInputData(); runAllTasks(); }
} 

void HistogramBase::buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
  completeGridObjectSetup(); actionsThatSelectTasks.push_back( getLabel() );
  if( !one_kernel_at_a_time ) { 
      unsigned hind = getNumberOfDerivatives();
      if( heights_index==2 && getPntrToArgument(arg_ends[hind])->getRank()>0 ) {
          tflags.assign(tflags.size(),1);
          for(unsigned i=0;i<tflags.size();++i) {
              if( fabs(getPntrToArgument(arg_ends[hind])->get( getTaskCode(i) ))<epsilon ) tflags[i]=0;
          }
      } else { tflags.assign(tflags.size(),1); }
      norm = static_cast<double>( tflags.size() ); 
  } else {
    std::vector<double> args( getNumberOfDerivatives() );
    double height=1.0; if( heights_index==2 ) height = getPntrToArgument(arg_ends[args.size()])->get();
    for(unsigned i=0; i<args.size(); ++i) args[i]=getPntrToArgument(i)->get();
    buildSingleKernel( tflags, height, args );
  }
}

void HistogramBase::performTask( const unsigned& current, MultiValue& myvals ) const {
  if( one_kernel_at_a_time ) {
    std::vector<double> args( getNumberOfDerivatives() ), der( getNumberOfDerivatives() );
    unsigned valout = getPntrToOutput(0)->getPositionInStream();
    gridobject.getGridPointCoordinates( current, args ); double vv = calculateValueOfSingleKernel( args, der );
    myvals.setValue( valout, vv );
    for(unsigned i=0; i<der.size(); ++i) { myvals.addDerivative( valout, i, der[i] ); myvals.updateIndex( valout, i ); }
  }
}

void HistogramBase::retrieveArgumentsAndHeight( const MultiValue& myvals, std::vector<double>& args, double& height ) const {
  std::vector<double> argsh( arg_ends.size()-1 ); retrieveArguments( myvals, argsh, 0 );
  height=1.0; if( heights_index==2 ) height = argsh[ argsh.size()-1 ];
  if( !unorm ) height = height / norm;
  for(unsigned i=0; i<args.size(); ++i) args[i]=argsh[i];
}

void HistogramBase::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                       const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( valindex==0 );
  if( one_kernel_at_a_time ) {
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
  if( one_kernel_at_a_time ) {
    if( getPntrToOutput(0)->forcesWereAdded() ) {
      unsigned valout = getPntrToOutput(0)->getPositionInStream(); double fforce = getPntrToOutput(0)->getForce( itask );
      for(unsigned i=0; i<getNumberOfDerivatives(); ++i) forces[i] += fforce*myvals.getDerivative( valout, i );
    }
    return;
  }
  std::vector<double> args( getNumberOfDerivatives() ); double height;
  retrieveArgumentsAndHeight( myvals, args, height );
  if( fabs(height)>epsilon ) {
    unsigned htask = 0; if( arg_ends[arg_ends.size()-1]-arg_ends[arg_ends.size()-2]>1 || getPntrToArgument(arg_ends[arg_ends.size()-2])->getRank()>0 ) htask=itask;
    if( getPntrToOutput(0)->forcesWereAdded() ) addKernelForces( heights_index, itask, args, htask, height, forces );
  }
}

}
}
