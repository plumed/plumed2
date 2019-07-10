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
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Communicator.h"

//+PLUMEDOC REWEIGHTING WHAM
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class Wham : 
public ActionWithValue,
public ActionWithArguments {
private:
  double thresh, simtemp;
  unsigned nreplicas;
  unsigned maxiter;
  bool wasUpdated;
  void calculateWeights();
public:
  static void registerKeywords(Keywords&);
  explicit Wham(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const { return 0; }
  void calculate(){}
  void apply(){}
  void update();
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(Wham,"WHAM")

void Wham::registerKeywords(Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.remove("ARG");
  keys.add("compulsory","ARG","the stored values for the bias");
  keys.add("compulsory","MAXITER","1000","maximum number of iterations for WHAM algorithm");
  keys.add("compulsory","WHAMTOL","1e-10","threshold for convergence of WHAM algorithm");
  keys.add("optional","TEMP","the system temperature.  This is not required if your MD code passes this quantity to PLUMED");
  keys.remove("NUMERICAL_DERIVATIVES");
}

Wham::Wham(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  wasUpdated(false)
{
  // Read in the temperature
  simtemp=0.; parse("TEMP",simtemp);
  if(simtemp>0) simtemp*=plumed.getAtoms().getKBoltzmann();
  else simtemp=plumed.getAtoms().getKbT();
  if(simtemp==0) error("The MD engine does not pass the temperature to plumed so you have to specify it using TEMP");
  // Now read in parameters of WHAM
  parse("MAXITER",maxiter); parse("WHAMTOL",thresh);
  if(comm.Get_rank()==0) nreplicas=multi_sim_comm.Get_size();
  comm.Bcast(nreplicas,0); addValue( getPntrToArgument(0)->getShape() ); 
  getPntrToOutput(0)->makeTimeSeries(); setNotPeriodic();
}

void Wham::calculateWeights() {
  // Retrieve the values that were stored for the biase
  std::vector<double> stored_biases( getPntrToArgument(0)->getNumberOfValues( getLabel() ) );
  for(unsigned i=0;i<stored_biases.size();++i) stored_biases[i] = getPntrToArgument(0)->get(i);
  // Get the minimum value of the bias
  double minv = *min_element(std::begin(stored_biases), std::end(stored_biases));
  // Resize final weights array
  plumed_assert( stored_biases.size()%nreplicas==0 );
  std::vector<double> final_weights( stored_biases.size() / nreplicas, 1.0 );
  if( getPntrToOutput(0)->getNumberOfValues( getLabel() )!=final_weights.size() ) {
      std::vector<unsigned> shape(1); shape[0]=final_weights.size(); getPntrToOutput(0)->setShape( shape );
  } 
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
    if( change<thresh ) {
        for(unsigned j=0; j<final_weights.size(); ++j) getPntrToOutput(0)->set( j, final_weights[j] );
        return; 
    }
  }
  error("Too many iterations in WHAM" );
}

void Wham::update() {
  if( getPntrToArgument(0)->getShape()[0]==0 ) return ;
  wasUpdated=true; calculateWeights();
}

void Wham::runFinalJobs() {
  if( !wasUpdated ) calculateWeights();
}

}
}
