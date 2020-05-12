/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
Calculate the weights for configurations using the weighted histogram analysis method.

Suppose that you have run multiple \f$N\f$ trajectories each of which was computed by integrating a different biased Hamiltonian. We can calculate the probability of observing
the set of configurations during the \f$N\f$ trajectories that we ran using the product of multinomial distributions shown below:
\f[
P( \vec{T} ) \propto \prod_{j=1}^M \prod_{k=1}^N (c_k w_{kj} p_j)^{t_{kj}}
\label{eqn:wham1}
\f]
In this expression the second product runs over the biases that were used when calculating the \f$N\f$ trajectories.  The first product then runs over the
\f$M\f$ bins in our histogram.  The \f$p_j\f$ variable that is inside this product is the quantity we wish to extract; namely, the unbiased probability of
having a set of CV values that lie within the range for the \f$j\f$th bin.

The quantity that we can easily extract from our simulations, \f$t_{kj}\f$, however, measures the number of frames from trajectory \f$k\f$ that are inside the \f$j\f$th bin.
To interpret this quantity we must consider the bias that acts on each of the replicas so the \f$w_{kj}\f$ term is introduced.  This quantity is calculated as:
\f[
w_{kj} =
\f]
and is essentially the factor that we have to multiply the unbiased probability of being in the bin by in order to get the probability that we would be inside this same bin in
the \f$k\f$th of our biased simulations.  Obviously, these \f$w_{kj}\f$ values depend on the value that the CVs take and also on the particular trajectory that we are investigating
all of which, remember, have different simulation biases.  Finally, \f$c_k\f$ is a free parameter that ensures that, for each \f$k\f$, the biased probability is normalized.

We can use the equation for the probability that was given above to find a set of values for \f$p_j\f$ that maximizes the likelihood of observing the trajectory.
This constrained optimization must be performed using a set of Lagrange multipliers, \f$\lambda_k\f$, that ensure that each of the biased probability distributions
are normalized, that is \f$\sum_j c_kw_{kj}p_j=1\f$.  Furthermore, as the problem is made easier if the quantity in the equation above is replaced by its logarithm
we actually chose to minimize
\f[
\mathcal{L}= \sum_{j=1}^M \sum_{k=1}^N t_{kj} \ln c_k  w_{kj} p_j + \sum_k\lambda_k \left( \sum_{j=1}^N c_k w_{kj} p_j - 1 \right)
\f]
After some manipulations the following (WHAM) equations emerge:
\f[
\begin{aligned}
p_j & \propto \frac{\sum_{k=1}^N t_{kj}}{\sum_k c_k w_{kj}} \\
c_k & =\frac{1}{\sum_{j=1}^M w_{kj} p_j}
\end{aligned}
\f]
which can be solved by computing the \f$p_j\f$ values using the first of the two equations above with an initial guess for the \f$c_k\f$ values and by then refining
these \f$p_j\f$ values using the \f$c_k\f$ values that are obtained by inserting the \f$p_j\f$ values obtained into the second of the two equations above.

Notice that only \f$\sum_k t_{kj}\f$, which is the total number of configurations from all the replicas that enter the \f$j\f$th bin, enters the WHAM equations above.
There is thus no need to record which replica generated each of the frames.  One can thus simply gather the trajectories from all the replicas together at the outset.
This observation is important as it is the basis of the binless formulation of WHAM that is implemented within PLUMED.

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
  bool buildsWeightStore() const override { return true; }
  void calculateWeights( const unsigned& nframes ) override;
  void clearData() override;
  double getLogWeight() override;
  double getWeight( const unsigned& iweight ) const override;
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
