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
#include "AverageBase.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "ActionRegister.h"

namespace PLMD {


void AverageBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys );
  ActionPilot::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.remove("ARG"); keys.use("UPDATE_FROM"); keys.use("UPDATE_UNTIL");
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the quantity being averaged");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
           "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.add("optional","LOGWEIGHTS","list of actions that calculates log weights that should be used to weight configurations when calculating averages");
}

AverageBase::AverageBase( const ActionOptions& ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  clearnextstep(false),
  firststep(true),
  clearnorm(false),
  data(getNumberOfArguments()),
  n_real_args(getNumberOfArguments())
{
  plumed_assert( keywords.exists("ARG") );
  std::vector<std::string> wwstr; parseVector("LOGWEIGHTS",wwstr);
  if( wwstr.size()>0 ) log.printf("  reweighting using weights from ");
  std::vector<Value*> arg( getArguments() ), biases; interpretArgumentList( wwstr, biases );
  for(unsigned i=0; i<biases.size(); ++i) {
    arg.push_back( biases[i] ); log.printf("%s ",biases[i]->getName().c_str() );
  }
  if( wwstr.size()>0 ) log.printf("\n");
  else log.printf("  weights are all equal to one\n");
  requestArguments( arg, false );

  // Read in clear instructions
  parse("CLEAR",clearstride);
  if( clearstride>0 ) {
    if( clearstride%getStride()!=0 ) error("CLEAR parameter must be a multiple of STRIDE");
    log.printf("  clearing average every %u steps \n",clearstride);
  }
}

void AverageBase::setupComponents( const unsigned& nreplicas ) { 
  unsigned nvals = getPntrToArgument(0)->getNumberOfValues( getLabel() );
  std::vector<unsigned> shape( 1 ); shape[0]=(clearstride / getStride() )*nvals*nreplicas; 
  for(unsigned j=0;j<n_real_args;++j) {
      if( getPntrToArgument(j)->getNumberOfValues( getLabel() )!=nvals ) error("all values input to store object must have same length");
      addComponent( getPntrToArgument(j)->getName(), shape ); 
      if( getPntrToArgument(j)->isPeriodic() ) { 
          std::string min, max; getPntrToArgument(j)->getDomain( min, max ); 
          componentIsPeriodic( getPntrToArgument(j)->getName(), min, max );
      } else componentIsNotPeriodic( getPntrToArgument(j)->getName() );
      getPntrToOutput(j)->makeTimeSeries();
  }
  // And create a component to store the weights
  addComponent( "logweights", shape ); componentIsNotPeriodic( "logweights" ); 
  getPntrToOutput( getNumberOfComponents()-1 )->makeTimeSeries();
}

std::string AverageBase::getStrideClearAndWeights() const {
  std::string stridestr; Tools::convert( getStride(), stridestr );
  std::string outstr = " STRIDE=" + stridestr;
  if( clearstride>0 ) {
      std::string clearstr; Tools::convert( clearstride, clearstr );
      outstr = " CLEAR=" + clearstr;
  }
  if( getNumberOfArguments()>n_real_args ) {
       outstr += " LOGWEIGHTS=" + getPntrToArgument(n_real_args)->getName();
       for(unsigned i=n_real_args+1; i<getNumberOfArguments(); ++i) outstr += "," + getPntrToArgument(i)->getName();
  }
  return outstr;
}

unsigned AverageBase::getNumberOfDerivatives() const {
  return getPntrToArgument(0)->getNumberOfDerivatives();
}

void AverageBase::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                                        std::vector<std::string>& max, std::vector<unsigned>& nbin,
                                        std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  plumed_dbg_assert( getNumberOfComponents()==1 && getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() );
  (getPntrToArgument(0)->getPntrToAction())->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, dumpcube );
}

void AverageBase::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  plumed_dbg_assert( getNumberOfComponents()==1 && getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() );
  (getPntrToArgument(0)->getPntrToAction())->getGridPointIndicesAndCoordinates( ind, indices, coords );
}

void AverageBase::getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const {
  plumed_dbg_assert( getNumberOfComponents()==1 && getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() );
  (getPntrToArgument(0)->getPntrToAction())->getGridPointAsCoordinate( ind, false, coords );
  if( coords.size()==(getPntrToOutput(0)->getRank()+1) ) coords[getPntrToOutput(0)->getRank()] = getPntrToOutput(0)->get(ind);
  else if( setlength ) {
    double val=getPntrToOutput(0)->get(ind);
    for(unsigned i=0; i<coords.size(); ++i) coords[i] = val*coords[i];
  }
}

void AverageBase::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

void AverageBase::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

void AverageBase::update() {
  // Resize values if they need resizing
  if( firststep ) { resizeValues(); firststep=false; }
  // Check if we need to accumulate
  if( (clearstride!=1 && getStep()==0) || !onStep() ) return;

  if( clearnextstep ) { 
      for(unsigned i=0;i<getNumberOfComponents();++i) {
          getPntrToOutput(i)->clearDerivatives(); getPntrToOutput(i)->set(0.0); 
      }
      if( clearnorm ) {
          for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToOutput(i)->setNorm(0.0);
      }
      clearnextstep=false; 
  }

  // Get the weight information
  double cweight=0.0; 
  if ( getNumberOfArguments()>n_real_args ) {
       for(unsigned i=n_real_args; i<getNumberOfArguments(); ++i) cweight+=getPntrToArgument(i)->get();
  }

  // Accumulate the data required for this round
  if( getPntrToArgument(0)->getRank()>0 && getPntrToArgument(0)->hasDerivatives() ) {
      accumulateGrid( cweight );
  } else {
      unsigned nvals = getPntrToArgument(0)->getNumberOfValues( getLabel() );
      for(unsigned i=0;i<nvals;++i) {
          for(unsigned j=0;j<n_real_args;++j) data[j] = getPntrToArgument(j)->get(i);
          accumulateValue( cweight, data );
      }
  }

  // Clear if required
  if( (clearstride>0 && getStep()%clearstride==0) ) clearnextstep=true;
}

void AverageBase::transferCollectedDataToValue( const std::vector<std::vector<double> >& mydata, const std::vector<double>& myweights ) {
  if( clearstride>0 ) return;
  std::vector<unsigned> shape(1); shape[0]=myweights.size();
  for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToOutput(i)->setShape( shape );

  for(unsigned i=0;i<myweights.size();++i) {
      for(unsigned j=0;j<getNumberOfComponents()-1;++j) { getPntrToOutput(j)->set( i, mydata[i][j] ); }
      getPntrToOutput(getNumberOfComponents()-1)->set( i, myweights[i] );
  }
}

}
