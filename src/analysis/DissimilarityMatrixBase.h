/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#ifndef __PLUMED_analysis_DissimilarityMatrixBase_h
#define __PLUMED_analysis_DissimilarityMatrixBase_h

#include "tools/Matrix.h"
#include "core/ActionPilot.h"

namespace PLMD {
namespace analysis {

class Analysis;

class DissimilarityMatrixBase : public ActionPilot {
private:
  Analysis* myanalysis_obj;
  unsigned nnodes;
  std::string fname, wfile;
  Matrix<double> mydissimilarities;
protected:
  void setDissimilarityMatrixElement( const unsigned& , const unsigned& , const double& );
public:
  static void registerKeywords( Keywords& keys );
  DissimilarityMatrixBase( const ActionOptions& ao );
  virtual double getDissimilarity( const unsigned& , const unsigned& ) const ;
  virtual void calculateDissimilarities(){ plumed_error(); }
  unsigned getNumberOfNodes() const ;
  void calculate(){};
  void apply(){};
  void update();
  void runFinalJobs();
};

inline
void DissimilarityMatrixBase::setDissimilarityMatrixElement( const unsigned& i, const unsigned& j, const double& v ){
  mydissimilarities(i,j) = v;
}

inline
double DissimilarityMatrixBase::getDissimilarity( const unsigned& i, const unsigned& j ) const {
  return mydissimilarities(i,j);
}

}
}
#endif
