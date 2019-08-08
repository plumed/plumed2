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
#include "core/AverageBase.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

//+PLUMEDOC REWEIGHTING ITRE
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class ITRE : 
public ActionWithValue,
public ActionWithArguments {
private:
  double thresh, simtemp;
  unsigned maxiter;
  bool wasUpdated;
  void calculateWeights();
public:
  static void registerKeywords(Keywords&);
  explicit ITRE(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const { return 0; }
  void calculate(){}
  void apply(){}
  void update();
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(ITRE,"ITRE")

void ITRE::registerKeywords(Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.remove("ARG");
  keys.add("compulsory","ARG","the stored values for the bias");
  keys.add("compulsory","MAXITER","1000","maximum number of iterations for WHAM algorithm");
  keys.add("compulsory","ITRETOL","1e-10","threshold for convergence of WHAM algorithm");
  keys.add("optional","TEMP","the system temperature.  This is not required if your MD code passes this quantity to PLUMED");
  keys.remove("NUMERICAL_DERIVATIVES");
}

ITRE::ITRE(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  wasUpdated(false)
{
  // Check argument is OK
  if( getNumberOfArguments()!=1 ) error("should only have one argument");
  Value* arg0 = getPntrToArgument(0);
  if( arg0->getRank()!=1 || arg0->hasDerivatives() || arg0->getName().find("logweights")==std::string::npos ) error("input should logweights from action that collects frames");
  // Setup the object that is collecting the data
  AverageBase* ab=dynamic_cast<AverageBase*>( arg0->getPntrToAction() );
  if( !ab ) error("input should be calculated by an average base object");
  ab->turnOnBiasHistory();
  // Read in the temperature
  simtemp=0.; parse("TEMP",simtemp);
  if(simtemp>0) simtemp*=plumed.getAtoms().getKBoltzmann();
  else simtemp=plumed.getAtoms().getKbT();
  if(simtemp==0) error("The MD engine does not pass the temperature to plumed so you have to specify it using TEMP");
  // Now read in parameters of WHAM
  parse("MAXITER",maxiter); parse("ITRETOL",thresh);
  // Now setup the value
  std::vector<unsigned> shape(1); shape[0] = getPntrToArgument(0)->getShape()[0];
  addValue( shape ); getPntrToOutput(0)->makeTimeSeries(); setNotPeriodic();
}

void ITRE::calculateWeights() {
  Value* arg0=getPntrToArgument(0); unsigned nvals=arg0->getShape()[0];
  // Set the shape of the output vector if we need to 
  if( getPntrToOutput(0)->getNumberOfValues( getLabel() )!=nvals ) {
      std::vector<unsigned> shape(1); shape[0]=nvals; getPntrToOutput(0)->setShape( shape );
  }
  Matrix<double> mymatrix( nvals, nvals );
  std::vector<double> in_weights( nvals ), new_weights( nvals ); 
  // Retrieve the weights from the input logweights
  in_weights[0] = mymatrix(0,0) = exp( arg0->get(0) );
  for(unsigned i=1;i<nvals;++i) {
      in_weights[i] = mymatrix(i,i) = exp( arg0->get( i*nvals + i ) );
      for(unsigned j=0;j<i;++j) mymatrix(i,j)=mymatrix(j,i) = exp( arg0->get( i*nvals + j ) );
  }
  // Now the iterative loop to calculate the ITRE weights
  for(unsigned iter=0; iter<maxiter; ++iter) {
    // This multiplies the vector in_weights by the matrix of weights -- the result is in new_weights
    mult( mymatrix, in_weights, new_weights ); 
 
    // Compute change in weights
    double change=0; 
    for(unsigned k=0; k<in_weights.size(); ++k) {
      double d = std::log( new_weights[k] / in_weights[k] ); change += d*d; in_weights[k] = new_weights[k];
    }
    change=epsilon;    // This is a tempory measure just to make this run and give something
    if( change<thresh ) {
        for(unsigned j=0; j<in_weights.size(); ++j) getPntrToOutput(0)->set( j, in_weights[j] );
        return; 
    }
  }
  error("Too many iterations in ITRE" );
}

void ITRE::update() {
  if( getPntrToArgument(0)->getShape()[0]==0 ) return ;
  wasUpdated=true; calculateWeights();
}

void ITRE::runFinalJobs() {
  if( !wasUpdated ) calculateWeights();
}

}
}
