/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#ifndef __PLUMED_analysis_Analysis_h
#define __PLUMED_analysis_Analysis_h

#include "core/ActionPilot.h"
#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "vesselbase/ActionWithVessel.h"

#define PLUMED_ANALYSIS_INIT(ao) Action(ao),Analysis(ao)

namespace PLMD {

class ReferenceConfiguration;

namespace analysis {

class AnalysisWithAnalysableOutput;

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new methods for analyzing the trajectory, within it there 
is \ref AddingAnAnalysis "information" as to how to go about implementing a new analysis method.

*/

class Analysis :
  public ActionPilot,
  public ActionAtomistic,
  public ActionWithArguments,
  public vesselbase::ActionWithVessel
  {
friend class AnalysisWithAnalysableOutput;
friend class ReadDissimilarityMatrix;
friend class LandmarkSelectionBase;
private:
/// Are we running only once for the whole trajectory
  bool single_run;
/// Are we treating each block of data separately
  bool nomemory;
/// Are we writing a checkpoint file
  bool write_chq;
/// Are we reusing data stored by another analysis action
  bool reusing_data;
/// Are we using data projected in a dimensionality reduction object
  bool dimred_data;
/// If we are reusing data are we ignoring the reweighting in that data
  bool ignore_reweight;
/// The Analysis action that we are reusing data from
  Analysis* mydatastash;
/// The frequency with which we are performing analysis
  unsigned freq;
/// Number of data point we need to run analysis
  unsigned ndata;
/// The temperature at which we are running the calculation
  double simtemp;
/// The temperature at which we want the histogram
  double rtemp;
/// Do we need the energy (are we reweighting at a different temperature)
  bool needeng;
/// The biases we are using in reweighting and the args we store them separately
  std::vector<Value*> biases;
/// The piece of data we are inserting
  unsigned idata;
/// The weights of all the data points
  std::vector<double> logweights;
/// Have we analyzed the data for the first time
  bool firstAnalysisDone;
/// The value of the old normalization constant
  double norm, old_norm;
/// Tempory vector to store values of arguments
  std::vector<double> current_args;
/// List of argument names 
  std::vector<std::string> argument_names;
/// The type of metric we are using to measure distances
  std::string metricname;
/// The checkpoint file
  OFile rfile;
/// Read in data from a file
  void readDataFromFile( const std::string& filename );
/// This retrieves the value of norm from the analysis action.
/// It is only used to transfer data from one analysis action to
/// another. You should never need to use it.  If you think you need it
/// you probably need getNormalization()
  double retrieveNorm() const ;
/// Get the metric if we are using malonobius distance and flexible hill
  std::vector<double> getMetric() const ;
/// The matrix containing the dissimiarities
  Matrix<double> mydissimilarities;
protected:
/// The dimensionality reduction object that we are reusing data from
  AnalysisWithAnalysableOutput* dimredstash;
/// Get a reference configuration 
  ReferenceConfiguration* getReferenceConfiguration( const unsigned& idata );
/// Get the squared dissimilarity between two reference configurations
  double getDissimilarity( const unsigned& i, const unsigned& j );
/// This is used to read in output file names for analysis methods.  When
/// this method is used and the calculation is not restarted old analysis
/// files are backed up.
  void parseOutputFile( const std::string& key, std::string& filename );
/// The data we are going to analyze
  std::vector<ReferenceConfiguration*> data;
/// Get the name of the metric we are using to measure distances
  std::string getMetricName() const ;
/// Return the number of arguments (this overwrites the one in ActionWithArguments)
  unsigned getNumberOfArguments() const;
/// Set the value of one fo the dissimilarities
  void setDissimilarityMatrixElement( const unsigned& , const unsigned& , const double& );
/// Retrieve the ith point
  void getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const ;
/// Returns true if argument i is periodic together with the domain 
  bool getPeriodicityInformation(const unsigned& i, std::string& dmin, std::string& dmax);
/// Return the normalization constant
  double getNormalization() const;
/// Return the set temperature (N.B. k_B T is what is returned by this function)
  double getTemp () const;
/// Are we analyzing each data block separately (if we are not this also returns the old normalization )
  bool usingMemory() const; 
/// Convert the stored log weights to proper weights
  void finalizeWeights( const bool& ignore_weights );
/// Overwrite ActionWithArguments getArguments() so that we don't return
/// the bias
  std::vector<Value*> getArguments();
/// Calculate the distance between stored snapshot iframe and stored snapshot jframe
  double getDistanceBetweenFrames( const unsigned& iframe, const unsigned& jframe, const bool& sq );
public:
  static void registerKeywords( Keywords& keys );
  Analysis(const ActionOptions&);
  ~Analysis();
/// Return the number of data points
  virtual unsigned getNumberOfDataPoints() const ;
/// Return the index of the data point in the base class
  virtual unsigned getDataPointIndexInBase( const unsigned& idata ) const ; 
/// Return the weight of the ith point
  virtual double getWeight( const unsigned& idata ) const ;
/// Calculate the dissimilarity between two configurations
  virtual void calcDissimilarity( const unsigned& , const unsigned& ){ plumed_error(); }
/// This returns the label of the object that contains the base data
  std::string getBaseDataLabel() const ;
  void prepare();
  void calculate();
  void update();
  void accumulate();
  virtual void performAnalysis()=0;
  void apply(){}
  void runFinalJobs();
  void runAnalysis();
  void lockRequests();
  void unlockRequests();
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ){ plumed_error(); }
  bool isPeriodic(){ plumed_error(); return false; }
  unsigned getNumberOfDerivatives(){ plumed_error(); return 0; }
  unsigned getRunFrequency() const ;
  bool runFinalAnalysisOnly() const ;
  virtual bool dissimilaritiesWereSet() const ;
};

inline
std::string Analysis::getMetricName() const {
  return metricname;
}

inline 
void Analysis::lockRequests(){
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
} 

inline
void Analysis::unlockRequests(){ 
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

inline
bool Analysis::usingMemory() const {
  if( !reusing_data ){
      return !nomemory;
  } else {
      return mydatastash->usingMemory();
  }
}

inline
unsigned Analysis::getNumberOfArguments() const {
  unsigned nargs=ActionWithArguments::getNumberOfArguments();
  return nargs - biases.size(); 
}

inline
double Analysis::retrieveNorm() const {
  return norm;
}

inline
std::vector<Value*> Analysis::getArguments(){
  std::vector<Value*> arg_vals( ActionWithArguments::getArguments() );
  for(unsigned i=0;i<biases.size();++i) arg_vals.erase(arg_vals.end()-1);
  return arg_vals;
}

inline
bool Analysis::runFinalAnalysisOnly() const {
  return single_run;
}

inline
void Analysis::setDissimilarityMatrixElement( const unsigned& i, const unsigned& j, const double& v ){
  plumed_dbg_assert( mydissimilarities(i,j)<epsilon );
  mydissimilarities(j,i)=mydissimilarities(i,j) = v;
}

}
}

#endif
