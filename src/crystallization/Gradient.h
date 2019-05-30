/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#ifndef __PLUMED_crystallization_Gradient_h
#define __PLUMED_crystallization_Gradient_h

#include "multicolvar/VolumeGradientBase.h"

namespace PLMD {
namespace crystallization {

class Gradient : public multicolvar::VolumeGradientBase {
  friend class GradientVessel;
private:
/// The value of sigma
  double sigma;
/// Number of quantities in use in this colvar
  unsigned vend, nquantities;
/// Number of bins in each direction
  std::vector<unsigned> nbins;
/// The type of kernel for the histogram
  std::string kerneltype;
public:
  static void registerKeywords( Keywords& keys );
  explicit Gradient(const ActionOptions&);
/// Get the number of quantities that are calculated each time
  virtual unsigned getNumberOfQuantities() const ;
/// Check on pbc - is it orthorhombic
  void setupRegions();
/// Calculate whats in the volume
  void calculateAllVolumes( const unsigned& curr, MultiValue& outvals ) const ;
};

inline
unsigned Gradient::getNumberOfQuantities() const {
  return nquantities;
}

}
}
#endif
