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
#ifndef __PLUMED_analysis_AnalysisWithDataCollection_h
#define __PLUMED_analysis_AnalysisWithDataCollection_h

#include "AnalysisBase.h"
#include "tools/PDB.h"

namespace PLMD {

class ReferenceConfiguration;

namespace analysis {

class ReadAnalysisFrames;

class AnalysisWithDataCollection : public AnalysisBase {
friend class EuclideanDissimilarityMatrix;
friend class ReadDissimilarityMatrix;
private:
/// Are we treating each block of data separately
  bool nomemory;
/// Are we writing a checkpoint file
  bool write_chq;
/// The temperature at which we are running the calculation
  double simtemp;
/// The temperature at which we want the histogram
  double rtemp;
/// The biases we are using in reweighting and the args we store them separately
  std::vector<Value*> biases;
/// The piece of data we are inserting
  unsigned idata;
/// The data we are going to analyze
  std::vector<ReferenceConfiguration*> data;
/// The weights of all the data points
  std::vector<double> logweights;
/// Have we analyzed the data for the first time
  bool firstAnalysisDone;
/// The value of the old normalization constant
  double norm, old_norm;
/// This allows you to output trajectory data with the projections
/// in EUCLIDEAN_DISSIMLAIRITIES and READ_DISSIMILARITY_MATRIX
  ReadAnalysisFrames* myframes;
/// Data is collected from the trajectory by passing it to this pdb.  These pdb
/// files are then read by the ReferenceConfigurations in data
  PDB mypdb;
/// The type of metric we are using to measure distances
  std::string metricname;
/// The checkpoint file --- really I would like to get rid of this and have some universal mechanism and a single file GT
  OFile rfile;
/// Perform the analysis -- we have a funciton as it is called from both runFinalJobs() and upate()
  void runAnalysis();
/// Read the checkpoint file  (this is used to read the nodes in readDissimilarityMatrix)
  void readCheckPointFile( const std::string& filename );
protected:
/// Return the temperature (used by Histogram)
  double getTemp() const { return simtemp; }
public:
  static void registerKeywords( Keywords& keys );
  AnalysisWithDataCollection(const ActionOptions&);
  ~AnalysisWithDataCollection();
/// Return the number of arguments (this overwrites the one in ActionWithArguments)
  unsigned getNumberOfArguments() const;
/// Return the number of data points
  virtual unsigned getNumberOfDataPoints() const ;
/// Return the index of the data point in the base class
  virtual unsigned getDataPointIndexInBase( const unsigned& idata ) const ;
/// Get the name of the metric we are using
  virtual std::string getMetricName() const ;
/// Return the weight of the ith point
  virtual double getWeight( const unsigned& idata ) const ;
/// Are we using memory in this calculation this affects the weights
  virtual bool usingMemory() const ;
/// Return the normalisation constant for the calculation
  virtual double getNormalization() const ;
/// By default dissimilarities are not set - they are only set in dissimilarity objects
  virtual bool dissimilaritiesWereSet() const ; 
/// This returns the label of the object that contains the base data
  std::string getBaseDataLabel() const ;
/// Get the ith data point
  virtual void getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const ;
/// Get a reference configuration (in dimensionality reduction this returns the projection)
  virtual ReferenceConfiguration* getReferenceConfiguration( const unsigned& idat, const bool& calcdist );
/// Now we can get the arguments from the right place
  const std::vector<Value*>    & getArguments() const ;
/// This ensures that the energy is stored if we are reweighting
  void prepare();
/// This stores the data and calls the analysis to be performed
  void update();
/// This does the analysis if it is to be done in the final step of the calculation
  void runFinalJobs();
};

inline
unsigned AnalysisWithDataCollection::getNumberOfDataPoints() const {
  if( !mydata ) return data.size();
  return AnalysisBase::getNumberOfDataPoints();
}

inline
std::string AnalysisWithDataCollection::getMetricName() const {
  if( !mydata ) return metricname; 
  return AnalysisBase::getMetricName();
}


inline
unsigned AnalysisWithDataCollection::getDataPointIndexInBase( const unsigned& idata ) const {
  if( !mydata ) return idata;
  return AnalysisBase::getDataPointIndexInBase( idata );
}

inline
bool AnalysisWithDataCollection::usingMemory() const {
  if( !mydata ) return !nomemory;
  return AnalysisBase::usingMemory();
}

inline
bool AnalysisWithDataCollection::dissimilaritiesWereSet() const { 
  if( !mydata ) return false; 
  return AnalysisBase::dissimilaritiesWereSet();
}

inline
double AnalysisWithDataCollection::getNormalization() const {
  if( !mydata ){
      if( nomemory || !firstAnalysisDone ) return norm;
      return ( 1. + norm/old_norm );
  } 
  return AnalysisBase::getNormalization();
}

inline
std::string AnalysisWithDataCollection::getBaseDataLabel() const {
  if( !mydata ) return getLabel();
  return AnalysisBase::getBaseDataLabel();
}

inline
unsigned AnalysisWithDataCollection::getNumberOfArguments() const {
  unsigned nargs=ActionWithArguments::getNumberOfArguments();
  return nargs - biases.size(); 
}

inline
const std::vector<Value*> & AnalysisWithDataCollection::getArguments() const {
  if( !mydata ) return ActionWithArguments::getArguments();
  return AnalysisBase::getArguments();
} 

}
}

#endif
