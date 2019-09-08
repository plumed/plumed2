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
#include "core/AverageBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace analysis {

class CollectFrames : public AverageBase {
private:
  unsigned ndata_for_norm, ndata;
  std::vector<double> data, off_diag_bias, allweights, posdata;
  std::vector<std::vector<double> > alldata;
public:
  static void registerKeywords( Keywords& keys );
  explicit CollectFrames( const ActionOptions& );
  void accumulateNorm( const double& cweight );
  void accumulateValue( const double& cweight, const std::vector<double>& dval );
  void accumulateAtoms( const double& cweight, const std::vector<Vector>& dir );
  unsigned getNumberOfStoredWeights() const ;
  void retrieveDataPoint( const unsigned& ipoint, const unsigned& jval, std::vector<double>& old_data );
  void storeRecomputedBias( const unsigned& ipoint, const unsigned& jframe, const double& data );
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(CollectFrames,"COLLECT_FRAMES")

void CollectFrames::registerKeywords( Keywords& keys ) {
  AverageBase::registerKeywords( keys ); ActionWithValue::useCustomisableComponents( keys );
  keys.add("optional","ARG","the data that you would like to collect to analyze later");
  keys.addOutputComponent("posx","ATOMS","these values store the x components of the atoms");
  keys.addOutputComponent("posy","ATOMS","these values store the y components of the atoms"); 
  keys.addOutputComponent("posz","ATOMS","these values store the z components of the atoms");
  keys.addOutputComponent("logweights","default","this value stores the logarithms of the weights of the stored configurations");
}

CollectFrames::CollectFrames( const ActionOptions& ao):
  Action(ao),
  AverageBase(ao),
  ndata_for_norm(0),
  ndata(0),
  data(n_real_args),
  posdata(3*getNumberOfAtoms())
{
  if( n_real_args>0 ) {
      if( getPntrToArgument(0)->hasDerivatives() && getPntrToArgument(0)->getRank()>0 ) {
          error("cannot collect grid input for later analysis -- if you need this email gareth.tribello@gmail.com");
      }
  }
  // Setup the components
  setupComponents( 1 ); 
}

void CollectFrames::accumulateNorm( const double& cweight ) {
  if( clearstride>0 ) {
      Value* bval=getPntrToOutput(getNumberOfComponents()-1);
      if( bval->getRank()==2 ) bval->set( ndata_for_norm*bval->getShape()[0] + ndata_for_norm, cweight );
      else bval->set( ndata_for_norm, cweight ); 
      ndata_for_norm++; if( getStep()%clearstride==0 ) ndata_for_norm=0; 
  } else allweights.push_back( cweight );  
}

void CollectFrames::accumulateValue( const double& cweight, const std::vector<double>& dval ) {
  // Now accumulate average
  if( clearstride>0 ) {
      for(unsigned j=0;j<dval.size();++j) getPntrToOutput(j)->set( ndata, dval[j] );
      ndata++; if( getStep()%clearstride==0 ) ndata=0;
      plumed_dbg_assert( ndata_for_norm==ndata );
  } else {
      alldata.push_back( dval );
  }
}

void CollectFrames::accumulateAtoms( const double& cweight, const std::vector<Vector>& dir ) {
  if( clearstride>0 ) {
      Vector thispos;
      for(unsigned i=0;i<dir.size();++i) {
          thispos = getReferencePosition(i) + dir[i];
          for(unsigned k=0;k<3;++k) getPntrToOutput(3*i+k)->set( ndata, thispos[k] );
      }
      ndata++; if( getStep()%clearstride==0 ) ndata=0; 
      plumed_dbg_assert( ndata_for_norm==ndata );
  } else {
      Vector thispos;
      for(unsigned i=0;i<dir.size();++i) {
          thispos = getReferencePosition(i) + dir[i];
          for(unsigned k=0;k<3;++k) posdata[3*i+k] = thispos[k];
      }
      alldata.push_back( posdata ); 
  }
}

unsigned CollectFrames::getNumberOfStoredWeights() const  {
  if( clearstride>0 ) return ndata_for_norm; 
  return allweights.size();
}

void CollectFrames::retrieveDataPoint( const unsigned& ipoint, const unsigned& jval, std::vector<double>& old_data ) {
  unsigned nvals=getPntrToArgument(0)->getNumberOfValues( getLabel() ); 
  if( clearstride>0 ) {
      plumed_dbg_assert( ipoint*nvals + jval < ndata );
      for(unsigned j=0;j<n_real_args;++j) old_data[j*nvals+jval] = getPntrToOutput(j)->get( ipoint*nvals + jval );
  } else {
      plumed_dbg_assert( ipoint*nvals + jval < alldata.size() );
      for(unsigned j=0;j<n_real_args;++j) old_data[j*nvals+jval] = alldata[ipoint*nvals + jval][j];

  }
}

void CollectFrames::storeRecomputedBias( const unsigned& ipoint, const unsigned& jframe, const double& biases ) {
  if( clearstride>0 ) {
      Value* bval=getPntrToOutput(getNumberOfComponents()-1);
      bval->set( bval->getShape()[0]*(ndata_for_norm+jframe) + ipoint + jframe, biases ); 
  } else {
      off_diag_bias.push_back( biases );
  }
}

void CollectFrames::runFinalJobs() {
  transferCollectedDataToValue( alldata, allweights, off_diag_bias );
}

}
}
