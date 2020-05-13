/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#ifndef __PLUMED_multicolvar_ActionVolume_h
#define __PLUMED_multicolvar_ActionVolume_h

#include "tools/HistogramBead.h"
#include "VolumeGradientBase.h"

namespace PLMD {
namespace multicolvar {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing a new way of defining a particular region of the simulation
box. You can use this to calculate the number of atoms inside that part or the average value of a quantity like the
coordination number inside that part of the cell.
*/

class ActionVolume : public VolumeGradientBase {
private:
/// Number of quantities in use in this colvar
  unsigned nquantities;
/// The value of sigma
  double sigma;
/// Are we interested in the area outside the colvar
  bool not_in;
/// The kernel type for this histogram
  std::string kerneltype;
protected:
  double getSigma() const ;
  std::string getKernelType() const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionVolume(const ActionOptions&);
/// Get the number of quantities that are calculated each time
  unsigned getNumberOfQuantities() const override;
/// Calculate whats in the volume
  void calculateAllVolumes( const unsigned& curr, MultiValue& outvals ) const override;
/// This calculates whether or not a particular is inside the box of interest
/// this is used for neighbour list with volumes
  bool inVolumeOfInterest( const unsigned& curr ) const ;
  virtual double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const=0;
  unsigned getCentralAtomElementIndex();
};

inline
unsigned ActionVolume::getNumberOfQuantities() const {
  return nquantities;
}

inline
double ActionVolume::getSigma() const {
  return sigma;
}

inline
std::string ActionVolume::getKernelType() const {
  return kerneltype;
}

inline
unsigned ActionVolume::getCentralAtomElementIndex() {
  return 1;
}

}
}
#endif
