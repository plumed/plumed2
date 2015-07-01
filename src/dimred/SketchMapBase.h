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
#ifndef __PLUMED_dimred_SketchMapBase_h
#define __PLUMED_dimred_SketchMapBase_h

#include "tools/SwitchingFunction.h"
#include "DimensionalityReductionBase.h"

namespace PLMD {
namespace dimred {

class SketchMapBase : public DimensionalityReductionBase {
private:
  bool reuse_hd, reuse_ld;
  SketchMapBase* smapbase;
  SwitchingFunction lowdf, highdf;
  std::vector<double> dtargets, ftargets;
protected:
  double mixparam;
public:
  static void registerKeywords( Keywords& keys );
  SketchMapBase( const ActionOptions& );
  void calculateProjections( const Matrix<double>& , Matrix<double>& );
  virtual void minimise( const Matrix<double>& , const Matrix<double>& , Matrix<double>& )=0;
  double transformLowDimensionalDistance( const double& val, double& df ) const ;
  double transformHighDimensionalDistance( const double& val, double& df ) const ;
  void setTargetDistance( const unsigned& idata, const double& dist );
  double calculateStress( const std::vector<double>& p, std::vector<double>& d );
};

inline
double SketchMapBase::transformLowDimensionalDistance( const double& val, double& df ) const {
  if( reuse_ld ) return smapbase->transformLowDimensionalDistance( val, df );
  double ans=1.0 - lowdf.calculate( val, df ); df*=-val; return ans;
}

inline
double SketchMapBase::transformHighDimensionalDistance( const double& val, double& df ) const {
  if( reuse_hd ) return smapbase->transformHighDimensionalDistance( val, df );
  double ans=1.0 - highdf.calculate( val, df ); df*=-val; return ans;
}

inline
void SketchMapBase::setTargetDistance( const unsigned& idata, const double& dist ){
  double df; dtargets[idata]=dist; ftargets[idata]=transformHighDimensionalDistance( dist, df );
}

}
}
#endif
