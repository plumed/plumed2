/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#ifndef __PLUMED_analysis_LandmarkSelectionBase_h
#define __PLUMED_analysis_LandmarkSelectionBase_h

#include "reference/ReferenceConfiguration.h"
#include "AnalysisWithLandmarks.h"

namespace PLMD {
namespace analysis {

class LandmarkSelectionOptions{
friend class LandmarkRegister;
friend class LandmarkSelectionBase;
private:
  std::vector<std::string> words;
  AnalysisWithLandmarks* action;
public:
  LandmarkSelectionOptions( const std::vector<std::string>& input, AnalysisWithLandmarks* myanalysis );
};

class LandmarkSelectionBase {
friend class AnalysisWithLandmarks;
friend class CopyAllFrames;
private:
/// Name of the method we are using for landmark selection
  std::string style;
/// The number of landmarks we are selecting
  unsigned nlandmarks;
/// The input to the landmark selection object
  std::vector<std::string> input;
/// A pointer to the underlying action
  AnalysisWithLandmarks* action;
/// How do we treat weights
  bool novoronoi, noweights;
protected:
/// Return the numbe of landmarks
  unsigned getNumberOfLandmarks() const ;
/// Return the communicator
  Communicator& getCommunicator();
/// Read a keywords from the input 
  template <class T>
  void parse(const std::string& ,T& );
/// Read a flag from the input
  void parseFlag(const std::string& key, bool& t);
/// Get the number of frames in the underlying action
  unsigned getNumberOfFrames() const;
/// Get the weight of the ith frame
  double getWeightOfFrame( const unsigned& );
/// Calculate the distance between the ith and jth frames
  double getDistanceBetweenFrames( const unsigned& , const unsigned&  );
/// Transfer frame i in the underlying action to the object we are going to analyze
  void selectFrame( const unsigned& , MultiReferenceBase* );
public:
  explicit LandmarkSelectionBase( const LandmarkSelectionOptions& lo );
  virtual ~LandmarkSelectionBase();
/// Check everything was read in
  void checkRead() const ;
/// Return a description of the landmark selection protocol
  std::string description();
/// Overwrite this to have a more descriptive output
  virtual std::string rest_of_description(){ return ""; };
/// Actually do landmark selection
  void selectLandmarks( MultiReferenceBase* );
  virtual void select( MultiReferenceBase* )=0;
};

inline
unsigned LandmarkSelectionBase::getNumberOfLandmarks() const {
  return nlandmarks;
}

inline
Communicator& LandmarkSelectionBase::getCommunicator(){
  return action->comm;
}

inline
unsigned LandmarkSelectionBase::getNumberOfFrames() const {
  return action->getNumberOfDataPoints();
}

template <class T>
void LandmarkSelectionBase::parse( const std::string& key, T& t ){
  bool found=Tools::parse(input,key,t);
  if(!found) plumed_merror("landmark seleciton style " + style + " requires " + key + " keyword");
}

}
}
#endif
