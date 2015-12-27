/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_gridtools_HistogramOnGrid_h
#define __PLUMED_gridtools_HistogramOnGrid_h

#include "GridVessel.h"

namespace PLMD {
namespace gridtools {

class HistogramOnGrid : public GridVessel {
private:
  double norm;
  bool store_normed;
  std::string kerneltype;
  std::vector<double> bandwidths;
  std::vector<unsigned> nneigh;
public:
  static void registerKeywords( Keywords& keys );
  explicit HistogramOnGrid( const vesselbase::VesselOptions& da );
  void setBounds( const std::vector<std::string>& smin, const std::vector<std::string>& smax,
                  const std::vector<unsigned>& nbins, const std::vector<double>& spacing );
  bool calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const ;
  void finish( const std::vector<double>& );
  bool applyForce(  std::vector<double>& forces ){ return false; }
  void addToNorm( const double& anorm );
  void setNorm( const double& snorm );
  double getNorm() const ;
  void switchOffNormalisation();
  void clear();
};

inline
void HistogramOnGrid::addToNorm( const double& anorm ){
  norm+=anorm;
}

inline
void HistogramOnGrid::setNorm( const double& snorm ){
  norm=snorm;
}

inline
double HistogramOnGrid::getNorm() const {
  return norm;
}

}
}
#endif
