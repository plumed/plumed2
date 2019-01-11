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
#ifndef __PLUMED_analysis_AverageVessel_h
#define __PLUMED_analysis_AverageVessel_h

#include "vesselbase/AveragingVessel.h"

namespace PLMD {
namespace analysis {

class AverageVessel : public vesselbase::AveragingVessel {
private:
  std::vector<double> domain;
public:
  /// keywords
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit AverageVessel( const vesselbase::VesselOptions& );
/// Set the size of the data vessel
  void resize();
/// This does nothing
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const ;
  std::string description() { return ""; }
/// Accumulate the average
  void accumulate( const double& weight, const double& val );
/// Get the average value
  double getAverage() const ;
};

}
}
#endif
