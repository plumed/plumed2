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

namespace PLMD {
namespace gridtools {

void HistogramBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG"); 
}

HistogramBase::HistogramBase(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  // Check all the values have the right size
  unsigned nvals=0; for(unsigned i=arg_ends[0];i<arg_ends[1];++i) nvals += getPntrToArgument(i)->getNumberOfValues();
  for(unsigned i=1;i<arg_ends.size()-1;++i) {
      unsigned tvals=0; for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) tvals += getPntrToArgument(j)->getNumberOfValues();
      if( nvals!=tvals ) error("mismatch between numbers of values in input arguments");
  }

  // Now construct the tasks
  if( actionInChain() ) {
      unsigned nrank = getPntrToArgument(0)->getShape()[0];
      for(unsigned i=1;i<getNumberOfArguments();++i) {
          if( arg_ends[i]!=i ) error("not sure if this sort of reshaping works");
          if( getPntrToArgument(i)->getShape()[0]!=nrank ) error("all arguments should have same shape");
      }
      for(unsigned i=0;i<nrank;++i) addTaskToList(i);
      if( distinct_arguments.size()>0 ){ 
          std::vector<std::string> alabels;
          for(unsigned i=0;i<getNumberOfArguments();++i){
              bool found=false; std::string mylab = (getPntrToArgument(i)->getPntrToAction())->getLabel();
              for(unsigned j=0;j<alabels.size();++j){
                  if( alabels[j]==mylab ){ found=true; break; }
              }
              if( !found ) alabels.push_back( mylab );
          }

          bool added=false;
          for(unsigned i=0;i<getNumberOfArguments();++i){
              // Add this function to jobs to do in recursive loop in previous action
              if( getPntrToArgument(i)->getRank()>0 ){
                  if( (getPntrToArgument(i)->getPntrToAction())->addActionToChain( alabels, this ) ){ added=true; break; }
              }
          }
          plumed_massert(added, "could not add action " + getLabel() + " to chain of any of its arguments");
      }
  } else {
      bool hasrank = getPntrToArgument(0)->getRank()>0;
      if( hasrank ) { 
          for(unsigned i=0;i<getNumberOfArguments();++i){ getPntrToArgument(i)->buildDataStore(); plumed_assert( getPntrToArgument(i)->getRank()>0 ); }
          for(unsigned i=0;i<nvals;++i) addTaskToList(i);
      } else {
          one_kernel_at_a_time=true; for(unsigned i=0;i<arg_ends.size();++i){ if( arg_ends[i]!=i ){ one_kernel_at_a_time=false; break; } }
          if( !one_kernel_at_a_time ) for(unsigned i=0;i<nvals;++i) addTaskToList(i);
      }
  }
}

void HistogramBase::addValueWithDerivatives( const std::vector<unsigned>& shape ){ 
  ActionWithValue::addValueWithDerivatives( shape ); setNotPeriodic();
  if( one_kernel_at_a_time ) {
      for(unsigned i=0;i<gridobject.getNumberOfPoints();++i) addTaskToList(i);
  }
}

unsigned HistogramBase::getNumberOfDerivatives() const {
  return arg_ends.size()-1;
}

void HistogramBase::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  gridobject.getGridPointCoordinates( ind, indices, coords );
}

void HistogramBase::calculate(){
  // Everything is done elsewhere
  if( actionInChain() ) return;
  // This is done if we are calculating a function of multiple cvs
  plumed_dbg_assert( getFullNumberOfTasks()>0 ); runAllTasks(); 
}

void HistogramBase::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  if( !actionInChain() && !one_kernel_at_a_time ) tflags.assign(tflags.size(),1);
  else {
     std::vector<double> args( arg_ends.size()-1 ); 
     for(unsigned i=0;i<arg_ends.size()-1;++i) args[i]=getPntrToArgument(arg_ends[i])->get();
     buildSingleKernel( tflags, args );
  }
}

void HistogramBase::performTask( const unsigned& current, MultiValue& myvals ) const {
  if( one_kernel_at_a_time ) { 
      std::vector<double> args( arg_ends.size()-1 ), der( arg_ends.size()-1 ); 
      unsigned valout = getPntrToOutput(0)->getPositionInStream(); 
      gridobject.getGridPointCoordinates( current, args ); double vv = calculateValueOfSingleKernel( args, der ); 
      myvals.setValue( valout, vv ); 
      for(unsigned i=0;i<der.size();++i){ myvals.addDerivative( valout, i, der[i] ); myvals.updateIndex( valout, i ); }
  }
}

void HistogramBase::gatherGridAccumulators( const unsigned& code, const MultiValue& myvals,
                                            const unsigned& bufstart, std::vector<double>& buffer ) const {
  if( one_kernel_at_a_time ) {
      unsigned istart = bufstart + (1+getNumberOfDerivatives())*code;
      unsigned valout = getPntrToOutput(0)->getPositionInStream(); buffer[istart] += myvals.get( valout );
      for(unsigned i=0;i<(arg_ends.size()-1);++i) buffer[istart+1+i] += myvals.getDerivative( valout, i );
      return;
  }
  std::vector<double> args( arg_ends.size()-1 );
  if( getPntrToArgument(0)->getRank()==2 ) {
      unsigned matind=getPntrToArgument(0)->getPositionInMatrixStash();
#ifdef DNDEBUG
      for(unsigned k=0;k<args.size();++k){
           unsigned amtind = getPntrToArgument(k)->getPositionInMatrixStash();
           plumed_dbg_assert( myvals.getNumberOfStashedMatrixElements(amtind)==myvals.getNumberOfStashedMatrixElements(matind) ); 
      }
#endif      
      for(unsigned j=0;j<myvals.getNumberOfStashedMatrixElements(matind);++j){
          unsigned jind = myvals.getStashedMatrixIndex(matind,j); 
          for(unsigned k=0;k<args.size();++k){
              unsigned amtind = getPntrToArgument(k)->getPositionInMatrixStash();
              plumed_assert( amtind==myvals.getStashedMatrixIndex(matind,j) );
              args[k] = myvals.getStashedMatrixElement( amtind, jind );
          }
          addKernelToGrid( args, bufstart, buffer );
      }
  } else { 
      retrieveArguments( myvals, args ); addKernelToGrid( args, bufstart, buffer );
  }
}

void HistogramBase::apply() { }

}
}
