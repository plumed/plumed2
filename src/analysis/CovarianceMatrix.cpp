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

//+PLUMEDOC REWEIGHTING COVARIANCE_MATRIX
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class CovarianceMatrix : 
public ActionWithValue,
public ActionWithArguments {
private:
  std::vector<double> weights;
  void buildTaskList();
  double getWeight( const unsigned& mycode ) const ;
public:
  static void registerKeywords(Keywords&);
  explicit CovarianceMatrix(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const override { return 0; }
  unsigned getNumberOfColumns() const override;
  void performTask( const unsigned& current, MultiValue& myvals ) const {}
  void buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags );
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals, const unsigned& bufstart, std::vector<double>& buffer ) const ;
  void calculate(){}
  void apply(){}
  void update();
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(CovarianceMatrix,"COVARIANCE_MATRIX")

void CovarianceMatrix::registerKeywords(Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","WEIGHTS","this keyword takes the label of an action that calculates a vector of values.  The elements of this vector "
           "are used as weights for the input data points.");
  keys.remove("NUMERICAL_DERIVATIVES");
}

CovarianceMatrix::CovarianceMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  unsigned ndata = arg_ends.size()-1;
  std::vector<std::string> weight_str; parseVector("WEIGHTS",weight_str);
  if( weight_str.size()>0 ) { 
    std::vector<Value*> weight_args, args( getArguments() ); interpretArgumentList( weight_str, weight_args );
    log.printf("  quantities used for weights are : %s ", weight_str[0].c_str() );
    for(unsigned i=1; i<weight_args.size(); ++i) log.printf(", %s", weight_str[i].c_str() );
    log.printf("\n");
  
    for(unsigned i=0; i<weight_args.size(); ++i) args.push_back( weight_args[i] );
    arg_ends.push_back( args.size() ); requestArguments( args, false );
  }
  buildTaskList(); std::vector<unsigned> shape(2); 
  shape[0]=shape[1]=ndata; addValue(shape); setNotPeriodic();
  getPntrToComponent(0)->alwaysStoreValues();
  // And convert eigenvector is a time series if it should be a time series  
  if( getPntrToArgument(0)->isTimeSeries() ) getPntrToComponent(0)->makeTimeSeries();
}

unsigned CovarianceMatrix::getNumberOfColumns() const {
  return arg_ends.size()-1;
}

void CovarianceMatrix::buildTaskList() {
  unsigned ndata=0;
  for(unsigned i=arg_ends[0];i<arg_ends[1];++i) ndata += getPntrToArgument(i)->getNumberOfValues();
  // Now check all data sets have same size
  for(unsigned i=0;i<arg_ends.size()-1;++i) {
      unsigned tval=0; for(unsigned j=arg_ends[i];j<arg_ends[i+1];++j) tval += getPntrToArgument(j)->getNumberOfValues();
      if( ndata>0 && tval!=ndata ) error("mismatch in sizes of input arguments for covariance matrix");
  } 
  // Now build the tasks
  for(unsigned i=0; i<ndata; ++i) addTaskToList(i);
}

double CovarianceMatrix::getWeight( const unsigned& mycode ) const {
  unsigned nt=0, nn=0, k=0, windex = arg_ends.size()-2;
  for(unsigned j=arg_ends[windex]; j<arg_ends[windex+1]; ++j) {
    unsigned nv = getPntrToArgument(j)->getNumberOfValues();
    nt += nv; if( mycode<nt ) { k=j; break; }
    nn += nv; k++;
  }
  return getPntrToArgument(k)->get( mycode - nn );
}

void CovarianceMatrix::buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
  actionsThatSelectTasks.push_back( getLabel() );
  double norm = 0; tflags.assign(tflags.size(),1);
  for(unsigned i=0;i<getFullNumberOfTasks();++i) {
      double ww = getWeight(i); norm += ww; if( std::fabs(ww)<epsilon ) tflags[i]=0;
  }
  if( weights.size()!=getFullNumberOfTasks() ) weights.resize( getFullNumberOfTasks() );
  // Setup the weights of all these points
  for(unsigned i=0;i<getFullNumberOfTasks();++i) weights[i] = getWeight(i) / norm; 
}

void CovarianceMatrix::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                          const unsigned& bufstart, std::vector<double>& buffer ) const {
   // Get the arguments
   plumed_dbg_assert( valindex==0 );
   unsigned nvals = arg_ends.size()-2; std::vector<double> argsh( arg_ends.size()-1 ); retrieveArguments( myvals, argsh, 0 ); 
   for(unsigned i=0;i<nvals;++i) {
       for(unsigned j=0;j<nvals;++j) buffer[bufstart + i*nvals + j] += weights[code]*argsh[i]*argsh[j]; 
   }
}

void CovarianceMatrix::update() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  runAllTasks();
}

void CovarianceMatrix::runFinalJobs() {
  if( skipUpdate() ) return;
  plumed_assert( !actionInChain() );
  // Need to create tasks here
  if( getFullNumberOfTasks()==0 ) { buildTaskList(); runAllTasks(); }
}

}
}
