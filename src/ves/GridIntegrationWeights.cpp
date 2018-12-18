/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "GridIntegrationWeights.h"

#include "tools/Grid.h"
#include "tools/File.h"
#include "tools/Exception.h"


namespace PLMD {
namespace ves {

void GridIntegrationWeights::getOneDimensionalIntegrationPointsAndWeights(std::vector<double>& points, std::vector<double>& weights, const unsigned int nbins, const double min, const double max, const std::string& weights_type) {
  double dx = (max-min)/(static_cast<double>(nbins)-1.0);
  points.resize(nbins);
  for(unsigned int i=0; i<nbins; i++) {points[i] = min + i*dx;}
  if(weights_type=="trapezoidal") {
    weights = getOneDimensionalTrapezoidalWeights(nbins,dx,false);
  }
  else {
    plumed_merror("getOneDimensionalIntegrationWeights: unknown weight type, the available type is trapezoidal");
  }
}


std::vector<double> GridIntegrationWeights::getOneDimensionalTrapezoidalWeights(const unsigned int nbins, const double dx, const bool periodic) {
  std::vector<double> weights_1d(nbins);
  for(unsigned int i=1; i<(nbins-1); i++) {
    weights_1d[i] = dx;
  }
  if(!periodic) {
    weights_1d[0]= 0.5*dx;
    weights_1d[(nbins-1)]= 0.5*dx;
  }
  else {
    // as for periodic arguments the first point should be counted twice as the
    // grid doesn't include its periodic copy
    weights_1d[0]= dx;
    weights_1d[(nbins-1)]= dx;
  }
  return weights_1d;
}


}
}
