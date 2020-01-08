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
#ifndef __PLUMED_analysis_ReadAnalysisFrames_h
#define __PLUMED_analysis_ReadAnalysisFrames_h

#include "AnalysisBase.h"
#include "bias/ReweightBase.h"

namespace PLMD {
namespace analysis {

class ReadAnalysisFrames : public AnalysisBase {
private:
/// The frequency with which to clear the data stash
  unsigned clearstride;
  bool clearonnextstep;
/// The list of argument names that we are storing
  std::vector<std::string> argument_names;
/// The list of atom numbers that we are storing
  std::vector<AtomNumber> atom_numbers;
/// The biases we are using in reweighting and the args we store them separately
  std::vector<Value*> weight_vals;
/// The object that calculates weights using WHAM
  bias::ReweightBase* wham_pointer;
/// The weights of all the data points
  bool weights_calculated;
  std::vector<double> logweights, weights;
/// The data that has been collected from the trajectory
  std::vector<DataCollectionObject> my_data_stash;
/// Calculate the weights of the various points from the logweights
  void calculateWeights();
public:
  static void registerKeywords( Keywords& keys );
  explicit ReadAnalysisFrames( const ActionOptions& ao );
  void update() override;
/// This does nothing
  void performAnalysis() override {}
/// This does nothing - it just ensures the final class is not abstract
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const override { plumed_error(); }
/// Get the number of data points
  unsigned getNumberOfDataPoints() const override;
/// Get the index of the data point
  unsigned getDataPointIndexInBase( const unsigned& idata ) const override;
/// Get the input arguments
  std::vector<Value*> getArgumentList() override;
/// Have dissimilarities between thses objects been calculated
  bool dissimilaritiesWereSet() const override;
/// How are dissimilarities calcualted is not known
  std::string getDissimilarityInstruction() const override;
/// Get the weight of one of the objects
  double getWeight( const unsigned& idat ) override;
/// Get the reference configuration
  DataCollectionObject & getStoredData( const unsigned& idata, const bool& calcdist ) override;
/// Get the list of atoms that are being stored
  const std::vector<AtomNumber>& getAtomIndexes() const override;
};

inline
unsigned ReadAnalysisFrames::getNumberOfDataPoints() const {
  return my_data_stash.size();
}

inline
unsigned ReadAnalysisFrames::getDataPointIndexInBase( const unsigned& idata ) const {
  return idata;
}

inline
bool ReadAnalysisFrames::dissimilaritiesWereSet() const {
  return false;
}

inline
double ReadAnalysisFrames::getWeight( const unsigned& idat ) {
  if( !weights_calculated ) calculateWeights();
  return weights[idat];
}

inline
DataCollectionObject & ReadAnalysisFrames::getStoredData( const unsigned& idata, const bool& calcdist ) {
  plumed_dbg_assert( idata<my_data_stash.size() );
  return my_data_stash[idata];
}

}
}
#endif
