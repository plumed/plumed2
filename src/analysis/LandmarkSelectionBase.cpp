/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "LandmarkSelectionBase.h"

namespace PLMD {
namespace analysis {

void LandmarkSelectionBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); ActionWithArguments::registerKeywords( keys );
  keys.use("ARG"); keys.add("compulsory","NLANDMARKS","the number of landmarks that you would like to select");
  ActionWithValue::useCustomisableComponents( keys );
}

LandmarkSelectionBase::LandmarkSelectionBase( const ActionOptions& ao ):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  jframe(0),
  nlandmarks(0),
  nvals(0)
{
  if( keywords.exists("NLANDMARKS") ) parse("NLANDMARKS",nlandmarks);
  log.printf("  selecting %u landmark points \n",nlandmarks);

  nvals = 0; 
  for(unsigned j=arg_ends[0];j<arg_ends[1];++j) {
      if( getPntrToArgument(j)->getRank()==1 ) nvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
      else if( getPntrToArgument(j)->getRank()==2 ) nvals += getPntrToArgument(j)->getShape()[1];
  }

  std::vector<unsigned> shape(1);
  for(unsigned i=0;i<arg_ends.size()-1;++i) { 
      unsigned tvals=0; 
      for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) {
          if( getPntrToArgument(j)->getRank()==1 ) tvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
          else if( getPntrToArgument(j)->getRank()==2 ) tvals += getPntrToArgument(j)->getShape()[1];
      }
      if( tvals!=nvals ) error("mismatch between sizes of input positions");
      // Add suitable argument
      Value* argi = getPntrToArgument(arg_ends[i]);
      if( argi->hasDerivatives() ) {
          error("cannot select landmarks for value " + argi->getName() );
      } else if( argi->getRank()==2 ) {
          shape.resize(2); shape[0] = nlandmarks; shape[1] = nvals;
      } else if( argi->getRank()==1 ) {
          shape.resize(1); shape[0] = nlandmarks;
      } else error("cannot select landmarks for value " + argi->getName() );
      addComponent( argi->getName(), shape );
      if( argi->isPeriodic() ) {
          std::string min, max; getPntrToArgument(arg_ends[i])->getDomain( min, max );
          componentIsPeriodic( getPntrToArgument(arg_ends[i])->getName(), min, max );
      } else componentIsNotPeriodic( getPntrToArgument(arg_ends[i])->getName() ); 
  }
}

void LandmarkSelectionBase::selectFrame( const unsigned& iframe ) {
  for(unsigned i=0;i<arg_ends.size()-1;++i) {
      Value* val0=getPntrToOutput(i);
      if( val0->getRank()==2 ) {
          Value* arg0 = getPntrToArgument(arg_ends[i]); 
          for(unsigned j=0;j<nvals;++j) val0->set( nvals*jframe + j, arg0->get(nvals*iframe+j) );
      } else if( val0->getRank()==1 ) val0->set( jframe, retrieveRequiredArgument( i, iframe ) );
  }
  jframe++; if( jframe==nlandmarks ) jframe=0;
}

void LandmarkSelectionBase::calculate() {
  if( skipCalculate() ) return;
  plumed_dbg_assert( !actionInChain() );
  selectLandmarks();
}

void LandmarkSelectionBase::update() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  selectLandmarks();
}

void LandmarkSelectionBase::runFinalJobs() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() ); nvals=0;
  for(unsigned j=arg_ends[0];j<arg_ends[1];++j) {
      if( getPntrToArgument(j)->getRank()==1 ) nvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
      else if( getPntrToArgument(j)->getRank()==2 ) nvals += getPntrToArgument(j)->getShape()[1];
  }
  std::vector<unsigned> shape(2); shape[0]=nlandmarks; shape[1]=nvals;
  for(unsigned i=0;i<arg_ends.size()-1;++i) {
      unsigned tvals=0;
      for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) {
          if( getPntrToArgument(j)->getRank()==1 ) tvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
          else if( getPntrToArgument(j)->getRank()==2 ) tvals += getPntrToArgument(j)->getShape()[1];
      }
      if( tvals!=nvals ) error("mismatch between sizes of input positions");
      Value* val0 = getPntrToComponent(i); if( val0->getRank()==2 ) { val0->setShape( shape ); } 
  }
  selectLandmarks();
}

}
}
