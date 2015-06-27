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
#ifndef __PLUMED_analysis_AnalysisWithAnalysableOutput_h
#define __PLUMED_analysis_AnalysisWithAnalysableOutput_h

#include "Analysis.h"

namespace PLMD {
namespace analysis {

class AnalysisWithAnalysableOutput : public Analysis {
private:
  std::vector<double> oweights;
protected:
  void setOutputWeights( const std::vector<double>& wwwin );
public:
  static void registerKeywords( Keywords& keys );
  AnalysisWithAnalysableOutput( const ActionOptions& );
  virtual unsigned getNumberOfOutputPoints() const=0;
  virtual ReferenceConfiguration* getOutputConfiguration( const unsigned& idata )=0;  
  virtual void getOutputForPoint( const unsigned& idata, std::vector<double>& point );
  double getOutputWeight( const unsigned& idata );
  virtual double getOutputDissimilarity( const unsigned& idata, const unsigned& jdata )=0;
};

inline
double AnalysisWithAnalysableOutput::getOutputWeight( const unsigned& idata ){
  plumed_dbg_assert( idata<getNumberOfOutputPoints() ); return oweights[idata]; 
}

}
}
#endif
