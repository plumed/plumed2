/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#ifndef __PLUMED_analysis_AnalysisBase_h
#define __PLUMED_analysis_AnalysisBase_h

#include "core/ActionPilot.h"
#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "vesselbase/ActionWithVessel.h"
#include "DataCollectionObject.h"

namespace PLMD {

class ReferenceConfiguration;

namespace analysis {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new methods for analyzing the trajectory. You can find
\ref AddingAnAnalysis "information" on how to use it to implement new analysis methods here.

*/

class AnalysisBase :
  public ActionPilot,
  public ActionWithValue,
  public ActionAtomistic,
  public ActionWithArguments,
  public vesselbase::ActionWithVessel
{
  friend class ReselectLandmarks;
  friend class ReadDissimilarityMatrix;
protected:
/// The Analysis action that we are reusing data from
  AnalysisBase* my_input_data;
public:
  static void registerKeywords( Keywords& keys );
  explicit AnalysisBase(const ActionOptions&);
/// These are required because we inherit from both ActionAtomistic and ActionWithArguments
  void lockRequests() override;
  void unlockRequests() override;
/// Return the number of data points
  virtual unsigned getNumberOfDataPoints() const ;
/// Return the index of the data point in the base class
  virtual unsigned getDataPointIndexInBase( const unsigned& idata ) const ;
/// Return the weight of the ith point
  virtual double getWeight( const unsigned& idata );
/// Get the name of the metric that is being used
  virtual std::string getMetricName() const ;
/// Are we using memory in this calculation this affects the weights of points
  virtual bool usingMemory() const ;
/// Return the normalisation constant for the calculation
  virtual double getNormalization() const ;
/// Ensures that dissimilarities were set somewhere
  virtual bool dissimilaritiesWereSet() const ;
/// Get the information on how dissimilarities were calculated for output PDB
  virtual std::string getDissimilarityInstruction() const ;
/// Get the squared dissimilarity between two reference configurations
  virtual double getDissimilarity( const unsigned& i, const unsigned& j );
/// Get the indices of the atoms that have been stored
  virtual const std::vector<AtomNumber>& getAtomIndexes() const ;
/// Overwrite getArguments so we get arguments from underlying class
  virtual std::vector<Value*> getArgumentList();
/// Get the list of argument names in the base
  std::vector<std::string> getArgumentNames();
/// Get a reference configuration (in dimensionality reduction this returns the projection)
  virtual DataCollectionObject& getStoredData( const unsigned& idata, const bool& calcdist );
/// This actually performs the analysis
  virtual void performAnalysis()=0;
/// These overwrite things from inherited classes (this is a bit of a fudge)
  bool isPeriodic() override { plumed_error(); return false; }
  unsigned getNumberOfDerivatives() override { plumed_error(); return 0; }
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ) override { plumed_error(); }
/// Calculate and apply do nothing all analysis is done during update step
  void calculate() override {}
  void apply() override {}
/// This will call the analysis to be performed
  void update() override;
/// This calls the analysis to be performed in the final step of the calculation
/// i.e. when use_all_data is true
  void runFinalJobs() override;
};

inline
void AnalysisBase::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

inline
void AnalysisBase::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

inline
unsigned AnalysisBase::getNumberOfDataPoints() const {
  return my_input_data->getNumberOfDataPoints();
}

inline
unsigned AnalysisBase::getDataPointIndexInBase( const unsigned& idata ) const {
  return my_input_data->getDataPointIndexInBase( idata );
}

inline
std::string AnalysisBase::getMetricName() const {
  return my_input_data->getMetricName();
}

inline
double AnalysisBase::getWeight( const unsigned& idata ) {
  return my_input_data->getWeight( idata );
}

inline
bool AnalysisBase::usingMemory() const {
  return my_input_data->usingMemory();
}

inline
double AnalysisBase::getNormalization() const {
  return my_input_data->getNormalization();
}

inline
bool AnalysisBase::dissimilaritiesWereSet() const {
  return my_input_data->dissimilaritiesWereSet();
}

inline
double AnalysisBase::getDissimilarity( const unsigned& i, const unsigned& j ) {
  return my_input_data->getDissimilarity( i, j );
}

inline
std::vector<Value*> AnalysisBase::getArgumentList() {
  return my_input_data->getArgumentList();
}

inline
DataCollectionObject& AnalysisBase::getStoredData( const unsigned& idata, const bool& calcdist ) {
  return my_input_data->getStoredData( idata, calcdist );
}

inline
const std::vector<AtomNumber>& AnalysisBase::getAtomIndexes() const {
  return my_input_data->getAtomIndexes();
}

inline
std::string AnalysisBase::getDissimilarityInstruction() const {
  return my_input_data->getDissimilarityInstruction();
}

}
}

#endif
