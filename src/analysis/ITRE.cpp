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
#include "core/CollectFrames.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <numeric>

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
  keys.add("compulsory","MAXITER","10","maximum number of iterations for iterative evaluation");
  keys.add("optional","TEMP","the system temperature");
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
  CollectFrames* ab=dynamic_cast<CollectFrames*>( arg0->getPntrToAction() );
  if( !ab ) error("input should be calculated by an average base object");
  ab->turnOnBiasHistory();
  // Read in the temperature
  simtemp=0.; parse("TEMP",simtemp);
  if(simtemp>0) simtemp*=plumed.getAtoms().getKBoltzmann();
  else simtemp=plumed.getAtoms().getKbT();
  if(simtemp==0) error("The MD engine does not pass the temperature to plumed so you have to specify it using TEMP");
  // Now read in parameters of WHAM
  parse("MAXITER",maxiter);
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
  Matrix<double> mymatrix( nvals, nvals ); // lower triangular exponential matrix
  std::vector<double> in_weights( nvals ); // weights
  // Retrieve the weights from the input logweights
  in_weights[0] = 1.0;
  for(unsigned i=0;i<nvals;++i) {
      in_weights[i] = 1.0;
      for(unsigned j=0;j<=i;++j) mymatrix(i,j) = exp( -simtemp*simtemp*arg0->get( i*nvals + j ) );
      for(unsigned j=i+1;j<nvals;++j) mymatrix(i,j)=0.0; // upper triangle empty
  }

  std::vector<double> denominator(nvals); // will store e[V(s(t),t) - c(t)]
  std::vector<double> normalization(nvals); // running integral of denominator
  // Now the iterative loop to calculate the ITRE weights
  for(unsigned iter=0; iter<maxiter; ++iter) {
    for (size_t index_1 = 0; index_1 < nvals; index_1++) {
      normalization[index_1] = 0.0; // I am doing this so that I am safe that normalization is 0
      denominator[index_1] = in_weights[index_1] / mymatrix(index_1,index_1) ;
    }
    std::partial_sum(denominator.begin(), denominator.end(), normalization.begin());
    // This multiplies the vector in_weights by the matrix of weights -- the result is in new_weights
    mult( mymatrix, denominator, in_weights );
    for (int index_1=0 ; index_1<nvals ; index_1++) in_weights[index_1] /= normalization[index_1];
  }
  for(unsigned j=0; j<in_weights.size(); ++j) getPntrToOutput(0)->set( j, -std::log(in_weights[j])/simtemp);
  return;
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
