/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2019 The plumed team
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
#include "ReweightBase.h"
#include "core/ActionRegister.h"
#include "tools/Communicator.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_WHAM
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightWham : public ReweightBase {
private:
  double thresh;
  unsigned nreplicas;
  unsigned maxiter;
  bool weightsCalculated;
  std::vector<double> stored_biases;
  std::vector<double> final_weights;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightWham(const ActionOptions&ao);
  bool buildsWeightStore() const { return true; }
  void calculateWeights( const unsigned& nframes );
  void clearData();
  double getLogWeight();
  double getWeight( const unsigned& iweight ) const ;
};

PLUMED_REGISTER_ACTION(ReweightWham,"REWEIGHT_WHAM")

void ReweightWham::registerKeywords(Keywords& keys ) {
  ReweightBase::registerKeywords( keys ); keys.remove("ARG");
  keys.add("compulsory","ARG","*.bias","the biases that must be taken into account when reweighting");
  keys.add("compulsory","MAXITER","1000","maximum number of iterations for WHAM algorithm");
  keys.add("compulsory","WHAMTOL","1e-10","threshold for convergence of WHAM algorithm");
}

ReweightWham::ReweightWham(const ActionOptions&ao):
  Action(ao),
  ReweightBase(ao),
  weightsCalculated(false)
{
  parse("MAXITER",maxiter); parse("WHAMTOL",thresh);
  if(comm.Get_rank()==0) nreplicas=multi_sim_comm.Get_size();
  comm.Bcast(nreplicas,0);
}

double ReweightWham::getLogWeight() {
  if( getStep()==0 ) return 1.0;  // This is here as first step is ignored in all analyses
  weightsCalculated=false;
  double bias=0.0; for(unsigned i=0; i<getNumberOfArguments(); ++i) bias+=getArgument(i);

  std::vector<double> biases(nreplicas,0.0);
  if(comm.Get_rank()==0) multi_sim_comm.Allgather(bias,biases);
  comm.Bcast(biases,0);
  for(unsigned i=0; i<biases.size(); i++) stored_biases.push_back( biases[i] );
  return 1.0;
}

void ReweightWham::clearData() {
  stored_biases.resize(0);
}

double ReweightWham::getWeight( const unsigned& iweight ) const {
  plumed_dbg_assert( weightsCalculated && iweight<final_weights.size() );
  return final_weights[iweight];
}

void ReweightWham::calculateWeights( const unsigned& nframes ) {
  if( stored_biases.size()!=nreplicas*nframes ) error("wrong number of weights stored");
  // Get the minimum value of the bias
  double minv = *min_element(std::begin(stored_biases), std::end(stored_biases));
  // Resize final weights array
  plumed_assert( stored_biases.size()%nreplicas==0 );
  final_weights.resize( stored_biases.size() / nreplicas, 1.0 );
  // Offset and exponential of the bias
  std::vector<double> expv( stored_biases.size() );
  for(unsigned i=0; i<expv.size(); ++i) expv[i] = exp( (-stored_biases[i]+minv) / simtemp );
  // Initialize Z
  std::vector<double> Z( nreplicas, 1.0 ), oldZ( nreplicas );
  // Now the iterative loop to calculate the WHAM weights
  for(unsigned iter=0; iter<maxiter; ++iter) {
    // Store Z
    for(unsigned j=0; j<Z.size(); ++j) oldZ[j]=Z[j];
    // Recompute weights
    double norm=0;
    for(unsigned j=0; j<final_weights.size(); ++j) {
      double ew=0;
      for(unsigned k=0; k<Z.size(); ++k) ew += expv[j*Z.size()+k]  / Z[k];
      final_weights[j] = 1.0 / ew; norm += final_weights[j];
    }
    // Normalize weights
    for(unsigned j=0; j<final_weights.size(); ++j) final_weights[j] /= norm;
    // Recompute Z
    for(unsigned j=0; j<Z.size(); ++j) Z[j] = 0.0;
    for(unsigned j=0; j<final_weights.size(); ++j) {
      for(unsigned k=0; k<Z.size(); ++k) Z[k] += final_weights[j]*expv[j*Z.size()+k];
    }
    // Normalize Z and compute change in Z
    double change=0; norm=0; for(unsigned k=0; k<Z.size(); ++k) norm+=Z[k];
    for(unsigned k=0; k<Z.size(); ++k) {
      Z[k] /= norm; double d = std::log( Z[k] / oldZ[k] ); change += d*d;
    }
    if( change<thresh ) { weightsCalculated=true; return; }
  }
  error("Too many iterations in WHAM" );
}

}
}
