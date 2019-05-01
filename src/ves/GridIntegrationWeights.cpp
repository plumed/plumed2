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

std::vector<double> GridIntegrationWeights::getIntegrationWeights(const Grid* grid_pntr, const std::string& fname_weights_grid, const std::string& weights_type) {
  std::vector<double> dx = grid_pntr->getDx();
  std::vector<bool> isPeriodic = grid_pntr->getIsPeriodic();
  std::vector<unsigned int> nbins = grid_pntr->getNbin();
  std::vector<std::vector<double> > weights_perdim;
  for(unsigned int k=0; k<grid_pntr->getDimension(); k++) {
    std::vector<double> weights_tmp;
    if(weights_type=="trapezoidal") {
      weights_tmp = getOneDimensionalTrapezoidalWeights(nbins[k],dx[k],isPeriodic[k]);
    }
    else {
      plumed_merror("getIntegrationWeights: unknown weight type, the available type is trapezoidal");
    }
    weights_perdim.push_back(weights_tmp);
  }

  std::vector<double> weights_vector(grid_pntr->getSize(),0.0);
  for(Grid::index_t l=0; l<grid_pntr->getSize(); l++) {
    std::vector<unsigned int> ind = grid_pntr->getIndices(l);
    double value = 1.0;
    for(unsigned int k=0; k<grid_pntr->getDimension(); k++) {
      value *= weights_perdim[k][ind[k]];
    }
    weights_vector[l] = value;
  }

  if(fname_weights_grid.size()>0) {
    Grid weights_grid = Grid(*grid_pntr);
    for(Grid::index_t l=0; l<weights_grid.getSize(); l++) {
      weights_grid.setValue(l,weights_vector[l]);
    }
    OFile ofile;
    ofile.enforceBackup();
    ofile.open(fname_weights_grid);
    weights_grid.writeToFile(ofile);
    ofile.close();
  }
  //
  return weights_vector;
}


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
