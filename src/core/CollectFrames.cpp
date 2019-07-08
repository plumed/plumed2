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

class CollectFrames : public AverageBase {
private:
  unsigned ndata;
  std::vector<double> data, allweights;
  std::vector<std::vector<double> > alldata;
public:
  static void registerKeywords( Keywords& keys );
  explicit CollectFrames( const ActionOptions& );
  void interpretDotStar( const std::string& ulab, unsigned& nargs, std::vector<Value*>& myvals );
  void accumulateData( const double& cweight );
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(CollectFrames,"COLLECT_FRAMES")

void CollectFrames::registerKeywords( Keywords& keys ) {
  AverageBase::registerKeywords( keys ); ActionWithValue::useCustomisableComponents( keys );
  keys.add("compulsory","ARG","the data that you would like to collect to analyze later");
  keys.addOutputComponent("weights","default","this value stores the weights of the stored configurations");
}

CollectFrames::CollectFrames( const ActionOptions& ao):
  Action(ao),
  AverageBase(ao),
  ndata(0),
  data(n_real_args)
{
  if( getPntrToArgument(0)->hasDerivatives() && getPntrToArgument(0)->getRank()>0 ) {
      error("cannot collect grid input for later analysis -- if you need this email gareth.tribello@gmail.com");
  } 
  unsigned nvals = getPntrToArgument(0)->getNumberOfValues( getLabel() );
  std::vector<unsigned> shape( 1 ); shape[0]=(clearstride / getStride() )*nvals; 
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
  addComponent( "weights", shape ); componentIsNotPeriodic( "weights" ); 
  getPntrToOutput( getNumberOfComponents()-1 )->makeTimeSeries();
}

void CollectFrames::interpretDotStar( const std::string& ulab, unsigned& nargs, std::vector<Value*>& myvals ) {
  for(unsigned i=0; i<getNumberOfComponents(); ++i) copyOutput(i)->interpretDataRequest( ulab, nargs, myvals, "" );
}

void CollectFrames::accumulateData( const double& cweight ) {
  unsigned nvals = getPntrToArgument(0)->getNumberOfValues( getLabel() ); 
  // Now accumulate average
  if( clearstride>0 ) {
      for(unsigned i=0;i<nvals;++i) {
          for(unsigned j=0;j<getNumberOfComponents()-1;++j) getPntrToOutput(j)->set( ndata, getPntrToArgument(j)->get(i) );
          getPntrToOutput(getNumberOfComponents()-1)->set( ndata, cweight ); ndata++;
      } 
      // Start filling the data set again from scratch
      if( getStep()%clearstride==0 ) { ndata=0; }
  } else {
      for(unsigned i=0;i<nvals;++i) {
          for(unsigned j=0;j<getNumberOfComponents()-1;++j) data[j] = getPntrToArgument(j)->get(i);
          allweights.push_back( cweight ); alldata.push_back( data );
      }
  }
}

void CollectFrames::runFinalJobs() {
  if( clearstride>0 ) return;
  std::vector<unsigned> shape(1); shape[0]=allweights.size(); 
  for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToOutput(i)->setShape( shape );

  for(unsigned i=0;i<allweights.size();++i) {
      for(unsigned j=0;j<getNumberOfComponents()-1;++j) { getPntrToOutput(j)->set( i, alldata[i][j] ); }
      getPntrToOutput(getNumberOfComponents()-1)->set( i, allweights[i] ); 
  }
}

}
