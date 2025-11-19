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
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"

//+PLUMEDOC REWEIGHTING WHAM
/*
Calculate the weights for configurations using the weighted histogram analysis method.

The example input below shows how this command is used.

```plumed
#SETTINGS NREPLICAS=4
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
rp: RESTRAINT ARG=phi KAPPA=50.0 AT=@replicas:{-3.00,-1.45,0.10,1.65}
# Get the bias on each of the replicas
rep: GATHER_REPLICAS ARG=rp.bias
# Merge the biases on each of the replicas into a single vector
all: CONCATENATE ARG=rep.*
# Collect all the bias values
col: COLLECT TYPE=vector ARG=all STRIDE=1
wham: WHAM ARG=col TEMP=300
DUMPVECTOR ARG=wham FILE=wham_data
```

As illustrated in the example above this command is used when you have run $N$ trajectories each of which was computed by integrating a different biased Hamiltonian. The input above calculates
the probability of observing the set of configurations during the $N$ trajectories that we ran using the product of multinomial distributions using the formula shown below:

$$
P( \underline{T} ) \propto \prod_{j=1}^M \prod_{k=1}^N (c_k w_{kj} p_j)^{t_{kj}}
$$

In this expression the second product runs over the biases that were used when calculating the $N$ trajectories.  The first product then runs over the
$M$ bins in our histogram.  The $p_j$ variable that is inside this product is the quantity we wish to extract; namely, the unbiased probability of
having a set of CV values that lie within the range for the $j$th bin.

The quantity that we can easily extract from our simulations, $t_{kj}$, however, measures the number of frames from trajectory $k$ that are inside the $j$th bin.
To interpret this quantity we must consider the bias that acts on each of the replicas so the $w_{kj}$ term is introduced.  This quantity is calculated as:

$$
w_{kj} = e^{+\beta V_k(x_j)}
$$

and is essentially the factor that we have to multiply the unbiased probability of being in the bin by in order to get the probability that we would be inside this same bin in
the $k$th of our biased simulations. Obviously, these $w_{kj}$ values depend on the value that the CVs take and also on the particular trajectory that we are investigating
all of which, remember, have different simulation biases.  Finally, $c_k$ is a free parameter that ensures that, for each $k$, the biased probability is normalized.

We can use the equation for the probability that was given above to find a set of values for $p_j$ that maximizes the likelihood of observing the trajectory.
This constrained optimization must be performed using a set of Lagrange multipliers, $\lambda_k$, that ensure that each of the biased probability distributions
are normalized, that is $\sum_j c_kw_{kj}p_j=1$.  Furthermore, as the problem is made easier if the quantity in the equation above is replaced by its logarithm
we actually chose to minimize

$$
L = \sum_{j=1}^M \sum_{k=1}^N t_{kj} \ln c_k  w_{kj} p_j + \sum_k\lambda_k \left( \sum_{j=1}^N c_k w_{kj} p_j - 1 \right)
$$

After some manipulations the following (WHAM) equations emerge:

$$
\begin{aligned}
p_j & \propto \frac{\sum_{k=1}^N t_{kj}}{\sum_k c_k w_{kj}} \\
c_k & =\frac{1}{\sum_{j=1}^M w_{kj} p_j}
\end{aligned}
$$

which can be solved by computing the $p_j$ values using the first of the two equations above with an initial guess for the $c_k$ values and by then refining
these $p_j$ values using the $c_k$ values that are obtained by inserting the $p_j$ values obtained into the second of the two equations above.

Notice that only $\sum_k t_{kj}$, which is the total number of configurations from all the replicas that enter the $j$th bin, enters the WHAM equations above.
There is thus no need to record which replica generated each of the frames.  One can thus simply gather the trajectories from all the replicas together at the outset.
This observation is important as it is the basis of the binless formulation of WHAM that is implemented within PLUMED.

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace wham {

class Wham :
  public ActionWithValue,
  public ActionWithArguments {
private:
  double thresh, simtemp;
  unsigned nreplicas;
  unsigned maxiter;
public:
  static void registerKeywords(Keywords&);
  explicit Wham(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
  void prepare() override ;
  void calculate() override ;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(Wham,"WHAM")

void Wham::registerKeywords(Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","scalar/vector/matrix","the stored values for the bias");
  keys.add("compulsory","MAXITER","1000","maximum number of iterations for WHAM algorithm");
  keys.add("compulsory","WHAMTOL","1e-10","threshold for convergence of WHAM algorithm");
  keys.add("optional","TEMP","the system temperature.  This is not required if your MD code passes this quantity to PLUMED");
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.setValueDescription("vector","the vector of WHAM weights to use for reweighting the elements in a time series");
}

Wham::Wham(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
  // Read in the temperature
  simtemp=getkBT();
  if(simtemp==0) {
    error("The MD engine does not pass the temperature to plumed so you have to specify it using TEMP");
  }
  // Now read in parameters of WHAM
  parse("MAXITER",maxiter);
  parse("WHAMTOL",thresh);
  if(comm.Get_rank()==0) {
    nreplicas=multi_sim_comm.Get_size();
  }
  comm.Bcast(nreplicas,0);
  addValue( getPntrToArgument(0)->getShape() );
  setNotPeriodic();
}

void Wham::prepare() {
  std::vector<std::size_t> shape(1);
  shape[0] = getPntrToArgument(0)->getNumberOfValues() / nreplicas;
  if( getPntrToComponent(0)->getNumberOfValues()!=shape[0] ) {
    getPntrToComponent(0)->setShape( shape );
  }
}

void Wham::calculate() {
  // Retrieve the values that were stored for the biase
  std::vector<double> stored_biases( getPntrToArgument(0)->getNumberOfValues() );
  for(unsigned i=0; i<stored_biases.size(); ++i) {
    stored_biases[i] = getPntrToArgument(0)->get(i);
  }
  // Get the minimum value of the bias
  double minv = *min_element(std::begin(stored_biases), std::end(stored_biases));
  // Resize final weights array
  plumed_assert( stored_biases.size()%nreplicas==0 );
  std::vector<double> final_weights( stored_biases.size() / nreplicas, 1.0 );
  // Offset and exponential of the bias
  std::vector<double> expv( stored_biases.size() );
  for(unsigned i=0; i<expv.size(); ++i) {
    expv[i] = exp( (-stored_biases[i]+minv) / simtemp );
  }
  // Initialize Z
  std::vector<double> Z( nreplicas, 1.0 ), oldZ( nreplicas );
  // Now the iterative loop to calculate the WHAM weights
  for(unsigned iter=0; iter<maxiter; ++iter) {
    // Store Z
    for(unsigned j=0; j<Z.size(); ++j) {
      oldZ[j]=Z[j];
    }
    // Recompute weights
    double norm=0;
    for(unsigned j=0; j<final_weights.size(); ++j) {
      double ew=0;
      for(unsigned k=0; k<Z.size(); ++k) {
        ew += expv[j*Z.size()+k]  / Z[k];
      }
      final_weights[j] = 1.0 / ew;
      norm += final_weights[j];
    }
    // Normalize weights
    for(unsigned j=0; j<final_weights.size(); ++j) {
      final_weights[j] /= norm;
    }
    // Recompute Z
    for(unsigned j=0; j<Z.size(); ++j) {
      Z[j] = 0.0;
    }
    for(unsigned j=0; j<final_weights.size(); ++j) {
      for(unsigned k=0; k<Z.size(); ++k) {
        Z[k] += final_weights[j]*expv[j*Z.size()+k];
      }
    }
    // Normalize Z and compute change in Z
    double change=0;
    norm=0;
    for(unsigned k=0; k<Z.size(); ++k) {
      norm+=Z[k];
    }
    for(unsigned k=0; k<Z.size(); ++k) {
      Z[k] /= norm;
      double d = std::log( Z[k] / oldZ[k] );
      change += d*d;
    }
    if( change<thresh ) {
      for(unsigned j=0; j<final_weights.size(); ++j) {
        getPntrToComponent(0)->set( j, final_weights[j] );
      }
      return;
    }
  }
  error("Too many iterations in WHAM" );
}

}
}
