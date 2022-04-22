/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#ifndef __PLUMED_gridtools_SphericalKDE_h
#define __PLUMED_gridtools_SphericalKDE_h

#include "HistogramBase.h"

namespace PLMD {
namespace gridtools {

class SphericalKDE : public HistogramBase {
private:
  double hh;
  std::vector<double> center;
  unsigned nbins;
  double von_misses_norm;
  double von_misses_concentration;
public:
  static void registerKeywords( Keywords& keys );
  explicit SphericalKDE(const ActionOptions&ao);
  void setupNeighborsVector() override;
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& out_nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void buildSingleKernel( const double& height, std::vector<double>& args );
  double calculateValueOfSingleKernel( const std::vector<double>& args, std::vector<double>& der ) const ;
  void addKernelToGrid( const double& height, const std::vector<double>& args, const unsigned& bufstart, std::vector<double>& buffer ) const ;
  void addKernelForces( const bool& height_has_derivative, const unsigned& itask, const std::vector<double>& args, const unsigned& htask, const double& height, std::vector<double>& forces ) const override;
};

}
}
#endif
