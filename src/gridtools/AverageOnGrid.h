/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#ifndef __PLUMED_gridtools_AverageOnGrid_h
#define __PLUMED_gridtools_AverageOnGrid_h

#include "HistogramOnGrid.h"

namespace PLMD {
namespace gridtools {

class AverageOnGrid : public HistogramOnGrid {
public:
  static void registerKeywords( Keywords& keys );
  explicit AverageOnGrid( const vesselbase::VesselOptions& da );
  void accumulate( const unsigned& ipoint, const double& weight, const double& dens, const std::vector<double>& der, std::vector<double>& buffer ) const override;
  void accumulateForce( const unsigned& ipoint, const double& weight, const std::vector<double>& der, std::vector<double>& intforce ) const override { plumed_error(); }
  double getGridElement( const unsigned& ipoint, const unsigned& jelement ) const override;
  unsigned getNumberOfComponents() const override;
  void getFinalForces( const std::vector<double>& buffer, std::vector<double>& finalForces ) override { plumed_error(); }
};

inline
unsigned AverageOnGrid::getNumberOfComponents() const {
  if( noderiv ) return nper - 1;
  return nper / ( dimension + 1 ) - 1;
}

}
}
#endif
